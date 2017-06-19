# (C) Copyright 2016 Isai Salas Gonzalez
# #
# #    This file is part of gdobap.
# #
# #    gdobap is free software: you can redistribute it and/or modify
# #    it under the terms of the GNU General Public License as published by
# #    the Free Software Foundation, either version 3 of the License, or
# #    (at your option) any later version.
# #
# #    gdobap is distributed in the hope that it will be useful,
# #    but WITHOUT ANY WARRANTY; without even the implied warranty of
# #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# #    GNU General Public License for more details.
# #
# #    You should have received a copy of the GNU General Public License
# #    along with gdobap.  If not, see <http://www.gnu.org/licenses/>.
#!/nas/longleaf/apps/perl/5.18.2/bin/perl
use strict;
use warnings;

my $file="Orthogroups.txt";
my $prefix=$ARGV[0];
my $outfile=join "",$prefix,"_",$file;
open (OUT,">",$outfile) || die "$outfile could not be opened\n";
open (IN,"<",$file) || die "$file could not be opened\n";
while (<IN>)
{
	chomp();
	my @ele=split " ",$_;
	my $mod=join "",$prefix,"_",$ele[0];
	print OUT $mod;
	for (my $i=1;$i <  scalar @ele;$i++)
	{
		print OUT " ",$ele[$i];
	}
	print OUT "\n";
}
close IN;
close OUT;


