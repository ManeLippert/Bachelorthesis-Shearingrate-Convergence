#!/bin/sh
# The script must be run from the directory containing the saved frames
# ($nameNNNNN.x)
# Convert the frames to 'fmt' format (fmt can be jpg, png)
fmt=jpeg
quality=50
# the 50 factor can vary between 40 and 60
fps=25
name=img
name=movie
# prefix for the movie file name
movieprefix=movie_out

# First I used APS version: x264 version but not playable on windows or quicktime
# Then I used Xvid but very large files
# Coudnt get any of the mpeg options with lavc to work

#
# read the image dimensions
# image width and height must be multiple of 16
#
for i in $name*.$fmt
do
  img=$i
  break
done

size=`identify $img | awk '{print $3}'`
width=${size/x*/}
height=${size/*x/}
height=${height/+*/}
#
# compute the optimal bitrate 
#       br = 50 * 25 * width * height / 256
#
# the 50 factor can vary between 40 and 60
#
obr=`expr $width \* $height \* $quality \* $fps / 256`

## clean temporary files that can interfere with the compression phase
rm -f divx2pass.log frameno.avi

##Resize
#width=800
#height=620

##See http://www.mplayerhq.hu/DOCS/HTML/en/menc-feat-selecting-codec.html
#codec="xvid"
codec="mpeg4"
#codec="h263"
#codec="h263p"
#codec="x264" ## Modern and small
#codec="msmpeg4"
#codec="msmpeg4v2"

#contfmt="avi"
moviesuffix="avi"
#moviesuffix="mp4"
#moviesuffix="mov"
#moviesuffix="flv"

## (MPEG 4)
## set the MPEG4 codec options - you may experiment!
opt="vbitrate=$obr:keyint=132:v4mv=yes:mbd=2:vqmin=3:lumi_mask=0.07:dark_mask=0.2:scplx_mask=0.1:tcplx_mask=0.1:naq:trell=yes"
mencoder -ovc lavc -lavcopts vcodec=$codec:vpass=1:$opt -mf type=$fmt:w=$width:h=$height:fps=$fps mf://\*.$fmt -nosound -o /dev/null
mencoder -ovc lavc -lavcopts vcodec=$codec:vpass=2:$opt -mf type=$fmt:w=$width:h=$height:fps=$fps mf://\*.$fmt -nosound -o $movieprefix.$moviesuffix
mencoder -ovc lavc -lavcopts vcodec=$codec:vpass=3:$opt -mf type=$fmt:w=$width:h=$height:fps=$fps mf://\*.$fmt -nosound -o $movieprefix.$moviesuffix

##XVID (MPEG 4 part 2 ASP)
#mencoder mf://\*.$fmt -of $contfmt -mf fps=$fps -ovc $codec -xvidencopts pass=1:bitrate=$obr -nosound -o $movieprefix.$moviesuffix
#mencoder mf://\*.$fmt -of $contfmt -mf fps=$fps -ovc $codec -xvidencopts pass=2:bitrate=$obr -nosound -o $movieprefix.$moviesuffix

##X264 (MPEG 4 part 10 -AVC)
## http://www.mplayerhq.hu/DOCS/HTML/en/menc-feat-x264.html
#opts="bitrate=5000:bframes=1:me=umh:partitions=all:trellis=1:qp_step=4:qcomp=0.7:direct_pred=auto:keyint=300"
#mencoder -of $contfmt -ovc $codec -x264encopts pass=1:$opts mf://\*.$fmt -nosound -o $movieprefix.$moviesuffix
#mencoder -of $contfmt -ovc $codec -x264encopts pass=2:$opts mf://\*.$fmt -nosound -o $movieprefix.$moviesuffix

## cleanup
rm -f divx2pass.log

# Container affects where it can be played
# For example quicktime cannot handle mpeg4 wrapped in avi.
# lacv will automatically do correct container for the suffix - for others, specify.
# if using x264 will need to be remuxed into a mp4 container.
## http://www.mplayerhq.hu/DOCS/HTML/en/menc-feat-quicktime-7.html

## I have found it impossible to make files made with mencoder which work in quicktime or WMP
## The formats are designed to be proprietrary and to those players are designed 
## to not work with open formats such as xvid.
## http://bemasc.net/wordpress/2010/02/02/no-you-cant-do-that-with-h264/
## http://www.rantroulette.com/2010/02/submarine-software-licenses-the-gotcha-of-video-codecs/