#eval "$(ssh-agent -s)" 
#ssh-add ~/.ssh/key_github

#!/bin/bash

# SSHエージェントがすでに起動しているか確認
if ! pgrep -u "$USER" ssh-agent > /dev/null; then
    echo "[log] Starting ssh-agent..."
    echo "[log] \"$(ssh-agent -s)\" "
    eval "$(ssh-agent -s)"
else
    echo "ssh-agent is already running."
fi

# SSHキーを追加 (必要ならパスを変更)
#SSH_KEY="$HOME/.ssh/key_github"
SSH_KEY="$HOME/.ssh/github_key2"
if [ -f "$SSH_KEY" ]; then
    echo "[log] ssh-add \"$SSH_KEY\" "
    ssh-add "$SSH_KEY"
    echo "SSH key added: $SSH_KEY"
else
    echo "SSH key not found: $SSH_KEY"
fi
