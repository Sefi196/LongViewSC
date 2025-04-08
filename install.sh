#!/bin/bash
# Step 1: Install R dependencies using devtools
echo "Installing ggtranscript R package..."
R -e "devtools::install_github('dzhang32/ggtranscript')"

# Step 1: Confirm successful installation
echo "Creating launcher script..."

# Define the launcher path
LAUNCHER_PATH="$HOME/bin/LongViewSC"

# Create the directory if it doesn't exist
mkdir -p "$(dirname "$LAUNCHER_PATH")"

# Create the launcher script
cat <<EOF > "$LAUNCHER_PATH"
#!/bin/bash
Rscript -e "shiny::runApp('$(pwd)/app.R', launch.browser=TRUE)"
EOF

# Make the launcher script executable
chmod +x "$LAUNCHER_PATH"

# Step 2: Add to PATH if it's not already included
if ! grep -q "$HOME/bin" "$HOME/.bashrc"; then
  echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bashrc"
  echo 'Added ~/bin to PATH. Please run the following command to apply changes:'
  echo 'source ~/.bashrc'
else
  echo '~/bin is already in your PATH.'
fi

# Final message
echo "You can now launch the app anytime by running: LongViewSC"
