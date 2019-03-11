#!/usr/bin/bash
#echo "test"が10回実行される。

for i in `seq 10` #1-10の配列
do
echo "test"
done

<< COMMENTOUT
echo "Hi, Jiro!"
echo "Hi, Saburo!"
COMMENTOUT


#!/bin/bash
#read -sp "Password: " pass
#tty -s && echo
#echo "<<$pass>>"

#!/bin/bash

HOST=対象ホスト名(IPアドレス)
USER=ログインユーザ
PASS=パスワード
TARGET_DIR=バックアップディレクトリ(ファイル)
BACKUP_DIR=保存先ディレクトリ

expect -c "
set timeout -1
spawn scp -Cpr ${TARGET_DIR} ${USER}@${HOST}:${BACKUP_DIR}
expect {
  \" Are you sure you want to continue connecting (yes/no)? \" {
    send \"yse\r\"
    expect \"password:\"
    send \"${PASS}\r\"
  } \"password:\" {
    send \"${PASS}\r\"
  }
}
interact
"
