$ErrorActionPreference = "Stop"

Write-Host "Checking for changes..."
$status = git status --porcelain

if ($status) {
    Write-Host "Changes detected. Staging all files..."
    git add .
    
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $commitMsg = "Auto-update: $timestamp"
    
    Write-Host "Committing with message: $commitMsg"
    git commit -m "$commitMsg"
    
    Write-Host "Pushing to remote..."
    git push
    
    Write-Host "Success! Changes pushed." -ForegroundColor Green
} else {
    Write-Host "No changes to push." -ForegroundColor Yellow
}
