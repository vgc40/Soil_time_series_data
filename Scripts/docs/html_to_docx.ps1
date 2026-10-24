# Convert the HTML documents written by write_storyboard_html.R to .docx
# (and a PDF copy for checking) with Microsoft Word.
#
# The HTML carries figure placeholders of the form [[IMG|<png path>|<width in>]];
# each is replaced by the embedded picture. Files are converted in a temporary
# folder outside OneDrive and then copied to Docs/.
#
# Run from the repository root:
#   powershell -ExecutionPolicy Bypass -File Scripts/docs/html_to_docx.ps1

$ErrorActionPreference = "Stop"
$docs = Join-Path (Get-Location) "Docs"
$work = Join-Path $env:TEMP "storyboard_docx"
New-Item -ItemType Directory -Force $work | Out-Null
$names = @("Storyboard_v4", "Analysis_notes_v4")

$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
$word.AutomationSecurity = 3   # msoAutomationSecurityForceDisable

try {
    foreach ($name in $names) {
        $html = Join-Path $work "$name.html"
        Copy-Item (Join-Path $docs "$name.html") $html -Force
        $docx = Join-Path $work "$name.docx"
        $pdf = Join-Path $work "$name.pdf"

        # Open(FileName, ConfirmConversions, ReadOnly, AddToRecentFiles)
        $doc = $word.Documents.Open($html, $false, $false, $false)
        $word.ActiveWindow.View.Type = 3       # wdPrintView

        $doc.PageSetup.PaperSize = 2           # wdPaperLetter
        foreach ($m in "TopMargin", "BottomMargin", "LeftMargin", "RightMargin") {
            $doc.PageSetup.$m = 72
        }

        # Replace figure placeholders with embedded pictures.
        while ($true) {
            $range = $doc.Content
            $find = $range.Find
            $find.ClearFormatting()
            $find.MatchWildcards = $true
            if (-not $find.Execute("\[\[IMG\|*\]\]")) { break }
            $parts = $range.Text.Trim("[", "]").Split("|")
            $path = $parts[1] -replace "/", "\"
            $width = [double]$parts[2] * 72
            $range.Text = ""
            $pic = $doc.InlineShapes.AddPicture($path, $false, $true, $range)
            $ratio = $pic.Height / $pic.Width
            $pic.Width = [Math]::Min($width, 468)
            $pic.Height = $pic.Width * $ratio
        }

        foreach ($table in @($doc.Tables)) {
            $table.AutoFitBehavior(2)          # wdAutoFitWindow
        }

        $doc.SaveAs2($docx, 16)                # wdFormatXMLDocument
        $doc.ExportAsFixedFormat($pdf, 17)     # wdExportFormatPDF
        $doc.Close($false)
        Copy-Item $docx, $pdf $docs -Force
        Write-Output "Wrote $(Join-Path $docs "$name.docx")"
    }
}
finally {
    $word.Quit()
    [System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}
