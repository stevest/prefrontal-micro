# !/bin/bash
kill -9 `ps -ef | grep stefanos | grep nrniv | grep -v grep | awk '{print$2}'`
