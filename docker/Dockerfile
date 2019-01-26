# Build from the jetscape-base image
FROM jdmulligan/jetscape-base:v1test

# Clone JETSCAPE and download external packages
RUN mkdir ./JETSCAPE && cd JETSCAPE \
&& git clone https://github.com/JETSCAPE/JETSCAPE.git . \
&& cd external_packages \
&& /bin/bash -c "source get_music.sh" \
&& /bin/bash -c "source get_lbtTab.sh" \
&& /bin/bash -c "source get_iSS.sh" \
&& /bin/bash -c "source get_freestream-milne.sh" \
&& cd .. && mkdir build \
&& cd build \
&& cmake .. -Dmusic=ON \
&& make -j8

CMD [“/bin/bash”]
