use LWP::Simple;
my $acc = @ARGV[0];
#@acc_array = split(/,/, $acc_list);
#print $acc_list;
#append [accn] field to each accession
#for ($i=0; $i < @acc_array; $i++) {
 #  $acc_array[$i] .= "[accn]";
#}
$acc .= "[accn]";
#join the accessions with OR
#$query = join('+OR+',@acc_array);

#assemble the esearch URL
$base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esearch.fcgi?db=nucleotide&term=$acc&usehistory=y";
#print $url, "\n";
#post the esearch URL
$output = get $url;
die "Couldn't get info" unless defined $output;

#print $output, "\n";
#parse WebEnv and QueryKey
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);

#assemble the efetch URL
$url = $base . "efetch.fcgi?db=nuccore&query_key=$key&WebEnv=$web";
$url .= "&rettype=uilist&retmode=text";
#print $url, "\n";
#post the efetch URL
$fasta = get($url);
die "Couldn't get info" unless defined $fasta;
print $fasta;