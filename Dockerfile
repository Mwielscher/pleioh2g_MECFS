# Uses Ubuntu 22.04 with CRAN binaries available via apt
FROM rocker/r2u:22.04

# System deps + Arrow C++ + R arrow binary + locales
RUN apt-get update && apt-get install -y --no-install-recommends \
    # R pkg binaries
    r-cran-arrow r-cran-data.table r-cran-dplyr r-cran-readr r-cran-fs r-cran-purrr r-cran-tibble r-cran-remotes \
    # build helpers (lightweight)
    git cmake pkg-config \
    # nice-to-have for HTTPS/Git
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    # locale tools
    locales \
 && rm -rf /var/lib/apt/lists/*

# Fix LC_* warning
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8

# Install pleioh2g from GitHub (arrow already satisfied via apt)
RUN R -e 'remotes::install_github("yjzhao1004/pleioh2g", upgrade = "never")'

# Working directory inside the container
WORKDIR /work

# Default to an interactive R session (override with `bash` or `Rscript your.R`)
CMD ["R"]
