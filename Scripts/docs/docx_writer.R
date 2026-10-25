# =============================================================================
# Minimal .docx writer (WordprocessingML) used by write_storyboard_docx.R
# =============================================================================
#
# Builds Word documents directly from R, with no Word automation. Content is
# written with light inline markup: <b>, <i>, and <span class="caption|note">
# become run formatting; other tags (links) are dropped and their text kept.


xml_escape <- function(x) {
  x %>% str_replace_all("&", "&amp;") %>% str_replace_all("<", "&lt;") %>%
    str_replace_all(">", "&gt;") %>% str_replace_all('"', "&quot;")
}

# Split marked-up text into Word runs.
runs <- function(text, size = NULL) {
  tokens <- str_split(text, "(?=<)|(?<=>)")[[1]]
  bold <- FALSE; italic <- FALSE; style <- "plain"
  out <- character()
  for (tok in tokens[tokens != ""]) {
    if (str_detect(tok, "^<b>$")) { bold <- TRUE; next }
    if (str_detect(tok, "^</b>$")) { bold <- FALSE; next }
    if (str_detect(tok, "^<i>$")) { italic <- TRUE; next }
    if (str_detect(tok, "^</i>$")) { italic <- FALSE; next }
    if (str_detect(tok, '^<span class="caption">$')) { style <- "caption"; next }
    if (str_detect(tok, '^<span class="note">$')) { style <- "note"; next }
    if (str_detect(tok, "^</span>$")) { style <- "plain"; next }
    if (str_detect(tok, "^</?[a-z][^>]*>$")) next
    props <- c(
      if (bold) "<w:b/>",
      if (italic) "<w:i/>",
      if (style == "note") '<w:color w:val="7F6000"/>',
      if (style == "caption" || !is.null(size)) sprintf('<w:sz w:val="%d"/>', if (is.null(size)) 19 else size)
    )
    out <- c(out, sprintf('<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r>',
                          if (length(props)) paste0("<w:rPr>", paste(props, collapse = ""), "</w:rPr>") else "",
                          xml_escape(tok)))
  }
  paste(out, collapse = "")
}

para <- function(text, style = NULL, extra_ppr = "", size = NULL) {
  ppr <- paste0(if (!is.null(style)) sprintf('<w:pStyle w:val="%s"/>', style), extra_ppr)
  sprintf("<w:p>%s%s</w:p>", if (nzchar(ppr)) paste0("<w:pPr>", ppr, "</w:pPr>") else "",
          runs(text, size))
}

h1 <- function(text) para(text, "Heading1")
h2 <- function(text) para(text, "Heading2")
h3 <- function(text) para(text, "Heading3")
p <- function(...) {
  text <- paste0(...)
  if (str_detect(text, 'class="note"')) {
    para(text, extra_ppr = '<w:shd w:val="clear" w:color="auto" w:fill="FFF5D6"/>')
  } else {
    para(text)
  }
}
b <- function(text) sprintf("<b>%s</b>", text)
i <- function(text) sprintf("<i>%s</i>", text)
ul <- function(items) {
  paste(map_chr(items, ~ para(.x, "ListBullet",
                              '<w:numPr><w:ilvl w:val="0"/><w:numId w:val="1"/></w:numPr>')),
        collapse = "")
}

# `rel` gives relative column widths (default equal); the table spans the
# 6.5-inch text width.
tbl <- function(df, rel = rep(1, ncol(df))) {
  widths <- floor(9360 * rel / sum(rel))
  cell <- function(text, width, header) {
    sprintf('<w:tc><w:tcPr><w:tcW w:w="%d" w:type="dxa"/>%s</w:tcPr>%s</w:tc>',
            width,
            if (header) '<w:shd w:val="clear" w:color="auto" w:fill="E7ECF5"/>' else "",
            para(if (header) b(text) else text, size = 18))
  }
  row <- function(values, header = FALSE) {
    paste0("<w:tr>", paste(map2_chr(as.character(values), widths, cell, header = header),
                           collapse = ""), "</w:tr>")
  }
  borders <- paste0("<w:tblBorders>",
                    paste(sprintf('<w:%s w:val="single" w:sz="4" w:space="0" w:color="999999"/>',
                                  c("top", "left", "bottom", "right", "insideH", "insideV")),
                          collapse = ""), "</w:tblBorders>")
  paste0(
    '<w:tbl><w:tblPr><w:tblW w:w="9360" w:type="dxa"/>', borders,
    '<w:tblLayout w:type="fixed"/><w:tblCellMar><w:left w:w="80" w:type="dxa"/><w:right w:w="80" w:type="dxa"/></w:tblCellMar></w:tblPr>',
    "<w:tblGrid>", paste(sprintf('<w:gridCol w:w="%d"/>', widths), collapse = ""), "</w:tblGrid>",
    row(names(df), header = TRUE),
    paste(apply(df, 1, row), collapse = ""),
    "</w:tbl>", para("")
  )
}

# Embedded figures: registered here, written to word/media by write_docx().
media <- new.env()
media$files <- character()

png_size <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  header <- readBin(con, "raw", 24)
  to_int <- function(bytes) sum(as.integer(bytes) * 256^(3:0))
  c(width = to_int(header[17:20]), height = to_int(header[21:24]))
}

img <- function(path, width_in = 6.5) {
  media$files <- c(media$files, path)
  k <- length(media$files)
  dims <- png_size(path)
  cx <- round(min(width_in, 6.5) * 914400)
  cy <- round(cx * dims[["height"]] / dims[["width"]])
  sprintf(paste0(
    '<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:drawing>',
    '<wp:inline distT="0" distB="0" distL="0" distR="0"><wp:extent cx="%d" cy="%d"/>',
    '<wp:docPr id="%d" name="Figure %d"/>',
    '<a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">',
    '<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">',
    '<pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">',
    '<pic:nvPicPr><pic:cNvPr id="%d" name="image%d.png"/><pic:cNvPicPr/></pic:nvPicPr>',
    '<pic:blipFill><a:blip r:embed="rIdImg%d"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>',
    '<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="%d" cy="%d"/></a:xfrm>',
    '<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic>',
    '</a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>'),
    cx, cy, k, k, k, k, k, cx, cy)
}

write_docx <- function(blocks, path) {
  tmp <- file.path(tempdir(), paste0("docx_", basename(path)))
  unlink(tmp, recursive = TRUE)
  dir.create(file.path(tmp, "_rels"), recursive = TRUE)
  dir.create(file.path(tmp, "word", "_rels"), recursive = TRUE)
  dir.create(file.path(tmp, "word", "media"), recursive = TRUE)
  ns <- paste0(
    'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" ',
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" ',
    'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"'
  )
  write <- function(text, ...) writeLines(text, file.path(tmp, ...), useBytes = TRUE, sep = "")
  decl <- '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'

  write(paste0(decl,
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
    '<Default Extension="xml" ContentType="application/xml"/>',
    '<Default Extension="png" ContentType="image/png"/>',
    '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>',
    '<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml"/>',
    '<Override PartName="/word/numbering.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.numbering+xml"/>',
    '</Types>'), "[Content_Types].xml")

  write(paste0(decl,
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>',
    '</Relationships>'), "_rels", ".rels")

  image_rels <- character()
  for (k in seq_along(media$files)) {
    file.copy(media$files[k], file.path(tmp, "word", "media", sprintf("image%d.png", k)))
    image_rels <- c(image_rels, sprintf(
      '<Relationship Id="rIdImg%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/image%d.png"/>', k, k))
  }
  write(paste0(decl,
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    '<Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>',
    '<Relationship Id="rIdNumbering" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/numbering" Target="numbering.xml"/>',
    paste(image_rels, collapse = ""), '</Relationships>'), "word", "_rels", "document.xml.rels")

  heading <- function(id, name, size, color, before) {
    sprintf(paste0('<w:style w:type="paragraph" w:styleId="%s"><w:name w:val="%s"/>',
                   '<w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:qFormat/>',
                   '<w:pPr><w:keepNext/><w:spacing w:before="%d" w:after="80"/><w:outlineLvl w:val="%d"/></w:pPr>',
                   '<w:rPr><w:b/><w:color w:val="%s"/><w:sz w:val="%d"/></w:rPr></w:style>'),
            id, name, before, as.integer(str_extract(id, "\\d")) - 1, color, size)
  }
  write(paste0(decl, '<w:styles ', ns, '>',
    '<w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:cs="Calibri" w:eastAsia="Calibri"/>',
    '<w:sz w:val="22"/></w:rPr></w:rPrDefault><w:pPrDefault><w:pPr><w:spacing w:after="120" w:line="264" w:lineRule="auto"/></w:pPr></w:pPrDefault></w:docDefaults>',
    '<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/></w:style>',
    heading("Heading1", "heading 1", 36, "1F3864", 240),
    heading("Heading2", "heading 2", 28, "1F3864", 360),
    heading("Heading3", "heading 3", 24, "2F5496", 240),
    '<w:style w:type="paragraph" w:styleId="ListBullet"><w:name w:val="List Bullet"/><w:basedOn w:val="Normal"/>',
    '<w:pPr><w:spacing w:after="60"/></w:pPr></w:style>',
    '</w:styles>'), "word", "styles.xml")

  write(paste0(decl, '<w:numbering ', ns, '>',
    '<w:abstractNum w:abstractNumId="0"><w:lvl w:ilvl="0"><w:start w:val="1"/><w:numFmt w:val="bullet"/>',
    '<w:lvlText w:val="•"/><w:lvlJc w:val="left"/><w:pPr><w:ind w:left="360" w:hanging="360"/></w:pPr></w:lvl></w:abstractNum>',
    '<w:num w:numId="1"><w:abstractNumId w:val="0"/></w:num></w:numbering>'), "word", "numbering.xml")

  write(paste0(decl, '<w:document ', ns, '><w:body>', paste(blocks, collapse = ""),
    '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/>',
    '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" w:header="720" w:footer="720" w:gutter="0"/>',
    '</w:sectPr></w:body></w:document>'), "word", "document.xml")

  files <- list.files(tmp, recursive = TRUE, all.files = TRUE)
  out <- normalizePath(path, winslash = "/", mustWork = FALSE)
  unlink(out)
  zip::zip(out, files = files, root = tmp)
  media$files <- character()
  invisible(out)
}
