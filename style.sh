#!/bin/bash
set -e

# Get script directory
# Source - https://stackoverflow.com/a/246128
# Posted by dogbane, modified by community. See post 'Timeline' for change history
# Retrieved 2026-05-20, License - CC BY-SA 4.0
SRC_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Check if astyle exists
# Source - https://stackoverflow.com/a/677212
# Posted by lhunath, modified by community. See post 'Timeline' for change history
# Retrieved 2026-05-20, License - CC BY-SA 4.0
if ! command -v astyle >/dev/null 2>&1
then
   echo "Artistic style could not be found."
   exit 1
fi

STYLE_OPTIONS="--suffix=none \
               --formatted \
               --indent=spaces=3 \
               --convert-tabs \
               --style=allman \
               --add-braces \
               --pad-header \
               --lineend=linux \
               --preserve-date \
               --max-continuation-indent=60 \
               --align-pointer=name \
               --max-code-length=80 \
               --break-after-logical"


FILES="$(find jabber/ tests/ apps/ -type f -name "*.cpp" -or -name "*.hpp")"

STYLE_COMMAND="astyle ${STYLE_OPTIONS} ${FILES}"

# Execute astyle.
OUTPUT="$($STYLE_COMMAND)"

EXIT_STATUS=0
# Throw error if anything was formatted.
if [ $(echo "$OUTPUT" | grep -c "Formatted") -gt 0 ];
then
   echo "Error: The following files have been styled:"
   echo "$(echo "$OUTPUT" | awk '{ print "       "$2}')"
   echo "Please commit these changes."
   EXIT_STATUS=1
else
   echo "Code is styled!"
fi
echo ""

# Check if line widths are all below 80.
LINE_WIDTH_INFO="$(echo "$FILES" | xargs wc -L --total=never)"
BAD_FILES=$(echo "$LINE_WIDTH_INFO" | awk '{ if ($1 > 80) print $0 }')
if [[ ! -z "$BAD_FILES" ]];
then
   echo "Error: The following files are > 80 chars wide:"
   echo "$BAD_FILES"
   echo "Please wrap at 80 chars, re-apply styling, and commit all changes."
   EXIT_STATUS=1
else
   echo "All files are <= 80 chars wide!"
fi

echo ""
if [[ $EXIT_STATUS -eq 0 ]];
then
   echo "Success!"
fi 

exit $EXIT_STATUS
