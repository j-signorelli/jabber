#!/bin/bash

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
               --indent=spaces=3 \
               --style=allman"

STYLE_COMMAND="astyle ${STYLE_OPTIONS} --recursive --dry-run \
                  ${SRC_DIR}/jabber/*.cpp,*.hpp \
                  ${SRC_DIR}/apps/*.cpp,*.hpp \
                  ${SRC_DIR}/tests/*.cpp,*.hpp"

# Execute astyle.
OUTPUT="$($STYLE_COMMAND)"
echo "$OUTPUT"
echo ""

# Throw error if anything was formatted.
if [ $(echo "$OUTPUT" | grep -c "Formatted") -gt 0 ];
then
   echo "Applied styling! Be sure to commit changes."
   exit 1
else
   echo "Code is styled."
   exit 0
fi
