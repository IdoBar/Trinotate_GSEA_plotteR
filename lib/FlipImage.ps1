$path = $args[0]

[Reflection.Assembly]::LoadWithPartialName("System.Windows.Forms"); 
$i = new-object System.Drawing.Bitmap $path

$i.RotateFlip("Rotate90FlipNone")

$i.Save($args[1],"png")