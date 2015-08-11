#!/bin/bash

ROOT=$(dirname $0)
cd $ROOT
OUTPUT="index.html" 

HTML_HEADER="<html><head><title>chemreac source distributions</title><body>"
HTML_FOOTER="</body></html>"

i=0
echo $HTML_HEADER > $OUTPUT
echo "Source distributions of chemreac (see <a href=\"github.com/chemreac/chemreac\">github.com/bjodah/chemreac</a>)." >> $OUTPUT
echo "<table><tr><th>Filename</th><th>Size</th><th>MD5</th></tr>" >> $OUTPUT
for i in `find "$ROOT" -maxdepth 1 -mindepth 1 -type f -iname "chemreac*.tar.gz"| sort -r`; do
    file=`basename "$i"`
    echo "    <tr><td><a href=\"./$file\">$file</a></td><td>$(du -h $file | cut -f1)</td><td>$(md5sum $file | head -c 32)</td></tr>" >> $OUTPUT
done
echo "</table>" >> $OUTPUT
#echo "this is a <a href=\"$(basename $0)\">statically generated</a> index." >> $OUTPUT
echo $HTML_FOOTER >> $OUTPUT
