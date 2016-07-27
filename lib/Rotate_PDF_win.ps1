$cpdf=dir -r ..\ -Filter cpdf.exe
& $cpdf.FullName -rotate 90 $args[0] AND -upright -o  $args[1]
echo "Saved rotated file into $args[1]"

