# !/bin/bash
kill -9 `ps -ef | grep $USER | grep nrniv | grep -v grep | awk '{print$2}'`
