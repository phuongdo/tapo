FROM pdoviet/tapo:latest

RUN  mv /home/sgeadmin/save/BioApps/tapo-v1.1.2 /home/sgeadmin/save/BioApps/tapo-v1.1.3 \
     && rm -r /home/sgeadmin/save/BioApps/tapo-v1.1.3/target \
     && cd /home/sgeadmin/save/BioApps/tapo-v1.1.3/ \
     && wget https://github.com/phuongdo/tapo/releases/download/v1.1.3-alpha.0/target.tar.gz \
     && tar -xzvf target.tar.gz

CMD [/bin/bash]
