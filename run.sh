#!/usr/bin/bash

input=$1
database=$2
if [[ -z $input ]] || [[ -z $database ]];
then
  echo "Please enter input fasta and database file."
  echo "bash $0 [input] [database]"
  exit
fi

if ! makeblastdb -h &> /dev/null
then
  echo "Blast+ could not be installed"
  exit
fi

function prefix() {
  local p=$1
  echo $(basename $p | perl -p -e 's/\.fa|fas|fasta|fna|txt$//')
}

src_p=$(dirname $0)
source $src_p/conf.sh

p_name=$(prefix $input)
s_name=$(prefix $database)
mkdir -p TEMP Result

cat $input | perl -n -e '
  s/\r?\n//;
  if (m/^>([^ ]+)/){
      $title = $1;
  }elsif (defined $title){
      $title_len{$title} += length($_);
  }
  END{
      for my $title (sort {$a cmp $b} keys %title_len){
          print "$title\t$title_len{$title}\n";
      }
  }
' >TEMP/$p_name.len

export len_f=TEMP/$p_name.len

makeblastdb -in $database -dbtype nucl -parse_seqids
blastn -query $input -db $database -outfmt 6 -task blastn > TEMP/$p_name.$s_name.out
cat TEMP/$p_name.$s_name.out |
  awk -v s_p=$similarity_p '$3 > s_p' |
  perl -nal -e '
    BEGIN{
      open my $f, "<", $ENV{len_f} or die $!;
      while (my $r = <$f>){
        chomp($r);
        my @t = split(/\s+/, $r);
        $len{$t[0]} = int($t[1]) if $r ne "";
      }
    }
    if ( exists $len{$F[0]} ){
      if ( ($F[7] - $F[6] + 1) / $len{$F[0]} > $ENV{length_p} ){
        print $_;
      }
    }
  ' |
  perl -nal -e '
    push @{ $hash{$F[1]}{$F[0]} }, [[$F[6], $F[7]], [$F[8], $F[9]]];
    END{
      for my $d (keys %hash){
        for my $q (keys %{$hash{$d}}){
          my $num = scalar(@{ $hash{$d}{$q} });
          if ( $keep_1 == 1 and $num > 1 ){
            my $m_v = 0;
            my $m_i = "F";
            for my $i (0..$num - 1){
              $a = $hash{$d}{$q}[$i];
              if ( $m_v < $a->[0][1] - $a->[0][0] ){
                $m_v = $a->[0][1] - $a->[0][0];
                $m_i = $i;
              }
            }
            $hash{$d}{$q} = [$hash{$d}{$q}[$m_i]];
          }
          my $v = $hash{$d}{$q}[0];
          print join("\t", $q, $d, $v->[0][0], $v->[0][1], $v->[1][0], $v->[1][1]);
        }
      }
    }
  '  |
  perl -nal -e '
    my $up_len = int($F[2]) - 1;
    my ($up_r, $up_l);
    if ( int($F[4]) < int($F[5]) ){
      $up_r = int($F[4]) - 1 - $up_len;
      $up_l = $up_r - $ENV{upstream} + 1;
    }else{
      $up_r = int($F[4]) + 1 + $up_len;
      $up_l = $up_r + $ENV{upstream} - 1;
    }
    unless ($up_l < 1){
      print join("\t", $F[1], join("-", "$up_l", "$up_r"));
    }
  ' >TEMP/$s_name.scale

export scale=TEMP/$s_name.scale
cat $database | perl -p -e 's/\r?\n//;s/^>(.+)$/>$1\n/;s/^>/\n>/' |
  perl -n -e '
    BEGIN{
      sub seq_comp_rev{return seq_com(seq_rev(shift))}
      sub seq_com{return shift =~ tr/AGTCagtc/TCAGtcag/r}
      sub seq_rev{ return reverse shift }
      open my $file_fh,"<",$ENV{scale} or die $!;
      while(<$file_fh>){
        s/\r?\n//;
        my @temp = split(/\s+/,$_);
        push @scales,[@temp];
      }
      close $file_fh;
    }
    s/\r?\n//;
    if (m/^>/){
      my $title = $_;
      my $sequence = <>;
      for my $line (@scales) {
        my $name = $line->[0];
        my @scale = @{$line}[1..scalar(@$line) - 1];
        if ($title =~ m/\b\Q$name\E\b/){
          for my $s (@scale){
            my @temp = split(/-/,$s);
            my ($start,$end) = ($temp[0],$temp[1]);
            my $len = abs($end - $start) + 1;
            my $seg;
            if ($start > $end){
                $seg = substr($sequence,$end - 1,$len);
                $seg = seq_comp_rev($seg);
            }else{
                $seg = substr($sequence,$start - 1,$len);
            }
            push @total , [">$name:$start-$end",$seg];
          }
        }
      }
    }
    END{
      for my $line (@total){
        printf "%s\n%s\n",$line->[0],$line->[1];
      }
    }
  ' > TEMP/$d_name.upstream.fa

cat TEMP/$d_name.upstream.fa | perl -nal -MList::Util -MData::Dumper -e '
  BEGIN{
    sub count {
      my $s = uc(shift);
      my %hash = ();
      map { $hash{$_}++ } split //, $s;
      my %return = map { $_ => $hash{$_} } qw/A T G C/;
      return \%return;
    }
    sub cosd {
      my ($l1, $l2) = @_;
      my $l = List::Util::min(scalar(@$l1), scalar(@$l2));
      my $cos = List::Util::sum( map {$l1->[$_] * $l2->[$_]} 0..$l - 1 ) / 
                (sqrt(List::Util::sum(map { int($l1->[$_]) ** 2 } 0..$l - 1)) * 
                 sqrt(List::Util::sum(map { int($l2->[$_]) ** 2 } 0..$l - 1)));
      return $cos;
    }
  }
  if (m/^>([^ ]+)/){
    my $title = $1;
    my $sequence = <>;
    chomp($hash{$title} = $sequence);
  }
  END{
    my %h = ();
    my %no_nr = ();
    for my $title (keys %hash){
      my $s = $hash{$title};
      if ( ! $h{ $s }++ ){
        my $c = count( $s );
        $no_nr{$title} = [[map {$c->{$_}} qw/A T C G/], $s];
      }
    }
    my @titles = keys %no_nr;
    my $length = scalar( @titles );
    my $matrix = [map {[map { 1 } 0..$length - 1]} 0..$length - 1];
    for my $i (0..$length - 2){
      for my $j ($i+1..$length - 1){
        my $i_d = $no_nr{$titles[$i]}->[0];
        my $j_d = $no_nr{$titles[$j]}->[0];
        $matrix->[$i][$j] = $matrix->[$j][$i] = cosd($i_d, $j_d);
      }
    }

    for my $i (0..$length - 1){
      my $p_num = scalar( grep { $_ < 0.98 } map { $matrix->[$i][$_] } @{$matrix->[$i]} );
      my $name = $titles[$i];
      if ( $p_num <= $ENV{isolate} and $length > 30 ){
        print STDERR ">$name\n$no_nr{$name}[1]";
      }else{
        print ">$name\n$no_nr{$name}[1]";
      }
    }
  }
' >Result/$s_name.upstream.nr.fa 2>Result/$s_name.upstream.discard.fa

rm -rf TEMP
