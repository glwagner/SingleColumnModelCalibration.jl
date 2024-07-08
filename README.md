# SingleColumnModelCalibration

Calibration of models for convective- and wind-driven mixing in the ocean surface boundary layer.

## How to run a calibration problem

To calibrate a closure in Oceananigans, take the following steps:

1. Download LES profile data. LES profiles can be found here: https://www.dropbox.com/scl/fo/gpxz01em3kmjcwuoafhnm/ABdBWTxkqSTbEycw9z2cbhc?rlkey=py5x1i7qvcmcjya0cue23k9tz&e=1&dl=0
This data must be copied into `data/profiles`. This can be done by copy/pasting the following lines from the root directory of the repo:

```bash
mkdir -p data
wget https://www.dropbox.com/scl/fi/8qwe9hf6wxc3w0crhzpml/profiles.zip\?rlkey\=sm8f7rhfokitzc7cgwinhip2k
mv profiles.zip\?rlkey=sm8f7rhfokitzc7cgwinhip2k data/profiles.zip
tar xvf data/profiles.zip
mv profiles data
```

