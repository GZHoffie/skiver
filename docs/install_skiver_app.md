# Installing the `skiver` app


## Ubuntu

Download the `AppImage` file from the [latest release of skiver](https://github.com/GZHoffie/skiver/releases/latest), either from the webpage or via `wget`.

```bash
# Replace this with the link to the latest release
wget https://github.com/GZHoffie/skiver/releases/download/v0.3.0/Skiver_0.3.0_amd64.AppImage
chmod a+x ./Skiver_0.3.0_amd64.AppImage

# Run with the following command
./Skiver_0.3.0_amd64.AppImage
```

### Troubleshoot

- **Wayland issue**: If `AppImage` fails to launch on a  system based on Wayland Display Server Protocol (instead of X11), try:

  ```
  LD_PRELOAD=/usr/lib/libwayland-client.so ./Skiver_0.3.0_amd64.appimage
  ```


## Windows

- Download `Skiver_<version>_x64-setup.exe` from the [latest release of skiver](https://github.com/GZHoffie/skiver/releases/latest), and open the file.
- A message will appear saying "Windows protected your PC". Click "More info" and "Run anyway".
- Follow the steps in the installer to complete the setup.