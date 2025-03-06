# Window-Grapher Sea-View-Factor

![izukyutokaido](https://github.com/user-attachments/assets/72ed6a3a-2d29-426c-858b-479f502b70cb)

This project demonstrates how to set up and use a Docker container for a Python pipeline to calculate sea-view-factor.

## Project Structure

```
my-docker-container
├── Dockerfile
├── docker-compose.yml
├── input
│   ├── N02-20_RailroadSection.geojson
         -> For Japan railway data, you can download the file from https://nlftp.mlit.go.jp/ksj/gml/datalist/KsjTmplt-N02-v2_3.html
│   ├── settings_calc.json
│   ├── settings_input.json
│   └── sample.tif
         -> For Japan DEM data, you can download the file from https://www.gsi.go.jp/kiban/
├── output
├── src
│   ├── main.py
│   ├── calc.py
│   ├── svf.py
│   ├── interpolate.py
│   ├── line.py
│   ├── raster.py
│   └── save.py
├── package.json
├── README.md
├── LICENSE.md
├── NOTICE.md
└── requirements.txt
```

## Getting Started

To get a copy of the project up and running on your local machine, follow these steps.

### Prerequisites

- Docker installed on your machine.

### Installation

1. Clone the repository:
   ```sh
   git clone https://github.com/window-grapher/sea-view-factor.git
   cd docker-svf
   ```

2. Build the Docker image:

   To use Docker Compose, follow these steps:

   ```sh
   docker-compose up --build
   ```

3. Start the services:

   You can skip this step if you have already done the step above.

   ```sh
   docker-compose up
   ```

   The `docker-compose.yml` file mounts the `/input` directory to `/app/input` and the `/output` directory to `/app/output` inside the container.

4. Access the running container:
   ```sh
   docker exec -it svf-container bash
   ```

   The `-it` flag is necessary to open the standard output.

### LICENSE

This source code is licensed by Apache 2.0. For detail, see the LICENSE.md and NOTICE.md.
FYI: https://infra.apache.org/licensing-howto.html
