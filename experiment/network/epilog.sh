echo "Epilog script"
echo "Calling matlab"
echo "First argument is \"${1}\""
echo "First argument is ${2}"
matlab -nojvm -nodisplay -nosplash -r "cd('${1}'); jobepilog('${2}');"; 
