#!/usr/bin/env bash

# This can be used as a RNALalifold post-processor:
# this script post-processes a set of Stockholm files by doing
#   - determine maximum covariation level
#   - run RNAz on the alignment
#   - create covariance model from mlocarna stk (optional)

rnaz=$(which RNAz)
alifoldz=$(which alifoldz.pl)
alifoldmaxcovar="$HOME/Perl/Scripts-Work/aln/alifold_maxcovar.pl"
removegaponly="$HOME/Perl/Scripts-Work/aln/remove-gaponly.pl"
reformataln="$HOME/Perl/Scripts-Work/aln/reformat_aln.pl"
stripaln="$HOME/Perl/Scripts-Work/aln/strip_aln.pl"
wd=$(pwd)
CM=OFF
STRIP=OFF
LOCARNATE=OFF
log="${wd}/pp_RNALalifold.log"
logcsv="${wd}/pp_RNALalifold.log.csv"
touch $log
touch $logcsv

for i in "$@"
do
  case $i in
    -l)          # make RNAz locarna-aware
      LOCARNATE=ON
      shift
      ;;
    -p=*)
      P="${i#*=}" # RNAz class probability
      shift
      ;;
    -s)
      STRIP=ON    # strip MSA, i.e. make non-redundant
      shift
      ;;
    -cm)
        CM=ON     # special aptamer treatment if OFF by default
        shift
        ;;
    *)
        # unknown option
        ;;
  esac
done

for w in $(ls *stk)
do
  wbn=$(basename $w .stk)
  echo "[I pp_RNALalifold] processing $wbn"
  # $removegaponly -a $w > $wbn.clean.aln 2> /dev/null
  #  $reformataln -a $w.clean.aln -i clustal -o clustal > $w.reformat.aln

  if [[ "$STRIP" == "ON" ]]; then
    $stripaln -a ${wbn}.stk -f S > ${wbn}.strip.stk
    cp -f ${wbn}.strip.stk ${wbn}.stk
  fi

  esl-reformat clustal ${wbn}.stk > ${wbn}.aln
  RNAZ_OPTIONS="-d"
  if [[ "$LOCARNATE" == "ON" ]]; then
    RNAZ_OPTIONS="${RNAZ_OPTIONS} -l"
  fi
  $rnaz ${RNAZ_OPTIONS} ${wbn}.aln > ${wbn}.rnaz.txt
  if [[ -s "${wbn}.rnaz.txt" ]]; then # grep SVM class probability
    rnazprob=$(grep probability ${wbn}.rnaz.txt | perl -ne 's/[a-zA-Z\-\:\s]+//g;print')
  else
    rnazprob="-1"
  fi

  $alifoldz -f -t 0.0 < ${wbn}.aln > ${wbn}.alifoldz.txt
  alifoldzscore=$(tail -n 1 ${wbn}.alifoldz.txt | grep -m 1 '.')

  RNAALIFOLD_OPTIONS="-t4 --aln --color -r --cfactor 0.6 --nfactor 0.5 -p --aln-EPS-cols=200 --aln-stk -f S"
  RNAalifold ${RNAALIFOLD_OPTIONS} < ${wbn}.stk > ${wbn}.alifold.out
  #mv RNAalifold_results.stk ${wbn}.RNAalifold_results.stk
  #mv alidot.ps ${wbn}.alidot.ps
  #mv alifold.out ${wbn}.alifold.out
  #mv alirna.ps ${wbn}.alirna.ps
  #mv aln.ps ${wbn}.aln.ps
  convert ${wbn}_aln.ps ${wbn}_aln.pdf
  convert ${wbn}_ss.ps ${wbn}_ss.pdf

  if [[ "$CM" == ON ]] && [[ $( echo "$rnazprob >= 0.9" | bc -l) -eq 1  ]]
  then
    echo "Building CM for ${wbn} (rnazprob=$rnazprob)"
    cmbuild ${wbn}.cm ${wbn}.stk > ${wbn}.cmbuild.out 2> ${wbn}.cmbuild.err
    cmcalibrate ${wbn}.cm > ${wbn}.cmcalibrate.out 2> ${wbn}.cmcalibrate.err
  fi

  #RNAplot -a ${wbn}.mlocarna.stk --aln --covar --aln-EPS-cols=300 -t 4 --auto-id --id-prefix ${wbn}.mlocarna


  $alifoldmaxcovar < ${wbn}_ali.out >> $log
  maxcovarval=$?
  echo "X $wbn maxcovar $maxcovarval $rnazprob $alifoldzscore" >> $log
  echo "$wbn;$maxcovarval;$rnazprob;$alifoldzscore" >> $logcsv
done
cd ..
