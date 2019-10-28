for x in A B C; do 
	echo Getting gene coordinates from $x
	awk '{if($3 == "gene"){gsub(";", " ", $9); print $1, $9, $4, $5}}' data/${x}.gff3 | sed 's/ID=//' | awk '{print $1, $2, $4, $5}' >  data/${x}.genes_positions.txt
done
