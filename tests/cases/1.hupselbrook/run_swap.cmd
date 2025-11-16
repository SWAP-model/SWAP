echo off
echo Running SWAP test case...
..\..\..\builddir_windows_native\swap.exe

echo SWAP run completed.
echo Press any key to continue... (automated)
echo. | set /p dummy=
rem Simulate pressing Enter for automated runs
echo.

rem Clean up temporary files
if exist swap.swp del swap.swp
if exist swap_swap.log del swap_swap.log  
if exist swap.ok del swap.ok
if exist result.* del result.*
if exist reruns.log del reruns.log

echo Temporary files cleaned up.
exit /b 0