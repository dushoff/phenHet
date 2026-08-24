use strict;
use 5.10.0;

while(<>){
	chomp;
	s/\$\$(.*\\begin\{align\})/$1/;
	s/(\\end\{align\}.*)\$\$/$1/;
	say;
}
