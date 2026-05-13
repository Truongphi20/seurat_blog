FROM satijalab/seurat:5.5.0

ENV DEBIAN_FRONTEND=noninteractive

# System deps
RUN apt-get update && apt-get install -y \
    curl git build-essential sudo \
    && rm -rf /var/lib/apt/lists/*

# Node (for mystmd)
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
    && apt-get install -y nodejs

RUN npm install -g mystmd@1.8.2

# rstudio user
ARG USERNAME=rstudio
RUN echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Install debug tools and dependencies
RUN apt-get update && apt-get install -y gdb cmake zlib1g-dev less

# Install bashdb for version 5.2
RUN apt-get update && apt-get install -y build-essential texinfo wget \
    && wget https://github.com/Trepan-Debuggers/bashdb/releases/download/5.2-1.2.0/bashdb-5.2-1.2.0.tar.gz \
    && tar -xzf bashdb-5.2-1.2.0.tar.gz \
    && cd bashdb-5.2-1.2.0 && ./configure && make && make install \
    && cd .. && rm -rf bashdb-5.2-1.2.0*

USER $USERNAME
WORKDIR /workspace

ENTRYPOINT ["bash", "-l"]