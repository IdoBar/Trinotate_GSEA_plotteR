dir -d *_output_plots -name | 
Foreach-Object { 
	$dirname = $_
	cd $dirname
	dir *_vert_27_07_16.pdf -name | 
		Foreach-Object { 
			$outfile=$_ + '_rotated.pdf'
			& ..\lib\cpdf.exe -rotate 90 $_ AND -upright -o  $outfile
			echo "Saved rotated file into $outfile"
		}
	cd ..
}
