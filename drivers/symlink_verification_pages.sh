#
# after generating the verification pages ("reports") symlink them into a
# group-specific subdirectory.  errors in the symlinking (e.g., if the symlink
# already exists) are sent to /dev/null.
#
groupdirs=`ls -d /home/luke/local/stringcheese/*group*`

for g in $groupdirs; do

  echo 'symlinking '$g' ...'

  bn=`basename $g`

  outdir=/home/luke/local/stringcheese/verification_pages/$bn

  [ ! -d $outdir ] && mkdir $outdir

  ln -s $g/*/*png $outdir/. 2>/dev/null

done

# print the output to stdout
for i in *group??_namenan; do echo $i: `ls $i/*png | wc -l` ; done
for i in *group???_namenan; do echo $i: `ls $i/*png | wc -l` ; done
for i in *group????_namenan; do echo $i: `ls $i/*png | wc -l` ; done
