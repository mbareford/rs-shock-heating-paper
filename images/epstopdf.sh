cd ../
find ./ -name "*.pdf" | grep -ve "/pdf/" | awk '{system("rm "$0)}'
find ./ -name "*.eps" | grep -ve "/pdf/" | awk '{system("epstool --copy --bbox "$0" --output "$0".b")}'
find ./ -name "*.eps.b" | grep -ve "/pdf/" | awk '{system("epstopdf "$0)}'
find . -name "*.eps.pdf" -exec sh -c 'mv -v "$0" "${0%.eps.pdf}.pdf"' '{}' \;
find . -type f -name "*.eps.b" -exec rm -f '{}' \;
