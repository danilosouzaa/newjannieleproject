#! /bin/sh

astyle --style=knf -s4 --indent=spaces=4 --indent-switches --indent-col1-comments \
        -z2 --convert-tabs --remove-brackets \
       *.c *.h



