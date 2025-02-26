. /home/once/general-venv/bin/activate
N="347_29190"
L="/home/once/Projects/Sites/static/tmp/$N.log"
cd /home/once/Projects/prime-gap
python misc/record_check.py --log-files logs/"$N"*.log | tee "$L"
date >> "$L"
