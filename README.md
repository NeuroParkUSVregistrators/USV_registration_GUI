# USV Registration GUI

A Julia-based graphical interface for **manual registration and inspection of mouse ultrasonic vocalizations (USVs)** from WaveSurfer `.h5` recordings.

The GUI provides interactive spectrogram inspection, manual USV bounding-box annotation, acoustic feature visualization, annotation saving, and automatic spectrogram-image export.

## Features

* Load WaveSurfer `.h5` recordings
* Display USV spectrograms
* Manually draw USV bounding boxes
* Move and resize existing annotations
* Delete selected annotations
* Save annotations as `.jld`
* Reopen and edit existing `.jld` annotations
* Automatically export individual USVs as JPEG images
* Display the current file and unsaved-change status
* Navigate recordings using keyboard and mouse controls
* Calculate acoustic features for individual and population-level USVs

Displayed acoustic features include:

* Duration
* Loudness
* Spectral purity
* Mean frequency
* Pitch variance

---

## Repository Structure

```text
USV_registration_GUI/
│
├── temp7.3.jl
├── usv.glade
│
├── Project.toml
├── Manifest.toml
│
└── Necessities/
    ├── acoustic.jl
    ├── GUItype.jl
    ├── gui_function_temp.jl
    ├── GDK_KEYmap.jl
    ├── wavesurfer.jl
    └── JPEGsaver_function_temp.jl
```

`temp7.3.jl` is the main GUI script.

---

## Requirements

The GUI is written in Julia and currently uses packages including:

```text
Gtk
GtkObservables
Graphics
Colors
Cairo
CairoMakie
JLD
HDF5
FFTW
DSP
LinearAlgebra
Distributions
Images
ImageTransformations
ImageMagick
FileIO
Interpolations
```

The current implementation is primarily developed and tested on Windows.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/NeuroParkUSVregistrators/USV_registration_GUI.git
cd USV_registration_GUI
```

Start Julia from the repository directory and activate the environment:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

If a `Manifest.toml` is provided, `Pkg.instantiate()` will install the package versions used by the project.

---

## Running the GUI

Run:

```julia
include("temp7.3.jl")
```

The GUI will open a file-selection dialog.

Select either:

```text
recording.h5
```

for a new registration, or

```text
recording.jld
```

to reopen a previously annotated recording.

The corresponding `.h5` recording must be available when loading a `.jld` annotation.

---

## Basic Workflow

```text
WaveSurfer recording (.h5)
          │
          ▼
     Load in GUI
          │
          ▼
 Inspect spectrogram
          │
          ▼
 Draw / edit USV boxes
          │
          ▼
        Save
       /    \
      ▼      ▼
   .jld     JPEGs
```

### 1. Load recording

Open a WaveSurfer `.h5` file.

The GUI calculates the spectrogram and displays the recording.

### 2. Register USVs

Draw bounding boxes around individual vocalizations.

### 3. Inspect annotations

The selected USV is shown separately together with its acoustic features.

### 4. Save

Use the **Save** button or:

```text
Ctrl + S
```

The annotation is saved as:

```text
recording.jld
```

and individual USV images are exported automatically.

---

## Controls

| Control            | Action                          |
| ------------------ | ------------------------------- |
| Left mouse drag    | Draw new USV                    |
| Modified left drag | Move selected USV               |
| Right mouse drag   | Resize selected USV             |
| Mouse wheel        | Zoom time scale                 |
| `←`                | Move backward through recording |
| `→`                | Move forward through recording  |
| `↑`                | Select next USV                 |
| `↓`                | Select previous USV             |
| `Ctrl + S`         | Save                            |
| `Ctrl + Z`         | Delete selected USV             |
| `Delete`           | Delete selected USV             |
| `Backspace`        | Delete selected USV             |

---

## File Status Indicator

The title bar shows the currently loaded filename.

For example:

```text
USV260801_001.h5
```

When annotations have been changed since the last save:

```text
USV260801_001.h5 *
```

The `*` indicates unsaved changes.

After saving, the `*` disappears.

---

## Annotation Files

Annotations are stored in JLD format.

Example:

```text
USV260801_001.jld
```

The annotation file stores USV bounding-box coordinates separately from Gtk-specific GUI objects, allowing annotations to be loaded again in later sessions.

The `.jld` filename corresponds to the original `.h5` recording:

```text
USV260801_001.h5
USV260801_001.jld
```

---

## JPEG Export

Saving also exports individual USV spectrograms.

Default output directory:

```text
USV260801_001_jpegs/
```

Default image size:

```text
128 × 128
```

The exported dataset can subsequently be used for image-based analyses such as dimensionality reduction, clustering, or neural-network models.

---

## Acoustic Measurements

For each manually registered USV, the GUI calculates:

| Measurement     | Description                                      |
| --------------- | ------------------------------------------------ |
| Duration        | Vocalization duration in ms                      |
| Loudness        | Mean signal power expressed in dB                |
| Spectral purity | Ratio describing concentration of spectral power |
| Mean frequency  | Mean dominant frequency in kHz                   |
| Pitch variance  | Variance of dominant frequency in log₁₀(Hz²)     |

Distributions of these measurements are displayed as histograms in the GUI.

---

## Current Version

### Temp7.3

Major changes:

* Added current filename to window title
* Added unsaved-change (`*`) indicator
* Improved file switching and GUI reset behavior
* Reorganized GUI state handling
* Changed dependencies to repository-relative paths
* Updated automatic JPEG export

See `Temp7.3_ReleaseNote` for details.

---

## Notes

The GUI currently assumes WaveSurfer-style HDF5 recordings and the acquisition structure expected by `wavesurfer.jl`.

Recording format, sampling rate, and channel configuration should therefore be checked before applying the GUI to data acquired with a different setup.

---

## Authors / Contributors

NeuroPark USV Registrators

---

## Repository

https://github.com/NeuroParkUSVregistrators/USV_registration_GUI
