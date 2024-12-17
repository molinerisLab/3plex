# Fail fast
set -o errexit
export PATH=/opt/mamba/bin:$PATH && \
export mamba_prefix=/opt/mamba
__conda_setup="$('/opt/mamba/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/mamba/etc/profile.d/conda.sh" ]; then
        . "/opt/mamba/etc/profile.d/conda.sh"
    else
        export PATH="/opt/mamba/bin:$PATH"
    fi
fi
mamba env update --quiet --file /root/3plex.yaml
