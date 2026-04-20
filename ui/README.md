# Skiver App
```
# From repo root
cargo build --release

# Then create the binaries dir and symlink
mkdir -p /home/zhenhao/skiver/ui/src-tauri/binaries
ln -sf /home/zhenhao/skiver/target/release/skiver \
/home/zhenhao/skiver/ui/src-tauri/binaries/skiver-x86_64-unknown-linux-gnu
```

```
# npm install
cd /home/zhenhao/skiver/ui
npm run tauri dev
```

Build the app with

```
cd /home/zhenhao/skiver/ui
npm run tauri build
```