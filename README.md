# USV Registration GUI — Temp8.1.2 Laptop

A laptop-optimized version of the Julia-based GUI for **manual and automatic registration of mouse ultrasonic vocalizations (USVs)** from WaveSurfer `.h5` recordings.

`temp8.1.2_laptop.jl` retains the automatic USV detection and manual correction workflow of Temp8.1, while reorganizing the GUI for use on **smaller laptop displays**.

---

## What is different from Temp8.1?

Temp8.1.2 Laptop primarily changes the **GUI layout and interaction**, rather than the USV detection workflow itself.

| Feature                            | Temp8.1 | Temp8.1.2 Laptop    |
| ---------------------------------- | ------- | ------------------- |
| Automatic USV detection            | ✓       | ✓                   |
| Manual box correction              | ✓       | ✓                   |
| `.jld` annotation saving           | ✓       | ✓                   |
| Automatic JPEG export              | ✓       | ✓                   |
| Standard GUI layout                | ✓       | —                   |
| Laptop-optimized layout            | —       | ✓                   |
| Responsive spectrogram canvas      | —       | ✓                   |
| Compact acoustic-feature plots     | —       | ✓                   |
| Spectrogram frequency-axis overlay | —       | ✓                   |
| Air-flow panel                     | ✓       | Removed             |
| Dedicated laptop `.glade` layout   | —       | `usv_vlaptop.glade` |

The laptop version therefore uses the same basic workflow:

```text
Load recording
      ↓
Automatic detection
      ↓
Manual inspection / correction
      ↓
Save annotations
```

but presents the information in a more compact interface.

---

# Installation

Clone the repository:

```bash
git clone https://github.com/NeuroParkUSVregistrators/USV_registration_GUI.git
cd USV_registration_GUI
git checkout Laptop_ver
```

Start Julia from the repository directory:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

The script loads its supporting functions from the `Necessities/` directory automatically.

---

# Running Temp8.1.2 Laptop

Run:

```julia
include("temp8.1.2_laptop.jl")
```

A file-selection dialog will appear.

Select either:

```text
recording.h5
```

to start a new registration, or:

```text
recording.jld
```

to reopen an existing annotation.

When loading a `.jld` file, the corresponding recording data must also be available.

---

# Basic Usage

## 1. Load a recording

Select a WaveSurfer `.h5` recording.

The recording is converted into a spectrogram and displayed in the main panel.

The current filename is shown in the window title.

If the annotation has been modified but not saved, an asterisk is added:

```text
recording.h5 *
```

---

## 2. Navigate through the recording

The visible time range can be controlled using the sliders or mouse wheel.

Keyboard controls are also available:

| Input       | Action                          |
| ----------- | ------------------------------- |
| `←`         | Move backward through recording |
| `→`         | Move forward through recording  |
| Mouse wheel | Change displayed time scale     |
| `↑`         | Select next USV                 |
| `↓`         | Select previous USV             |
| `Ctrl + S`  | Save                            |
| `Ctrl + Z`  | Delete selected USV             |
| `Delete`    | Delete selected USV             |
| `Backspace` | Delete selected USV             |

The laptop version handles these shortcuts at the main-window level so that interacting with sliders or other GUI elements does not unnecessarily capture keyboard focus.

---

# Automatic USV Detection

Press the **Detection** button to automatically identify candidate USVs.

The detection workflow is:

```text
Original spectrogram
        ↓
      Denoise
        ↓
Temporal USV detection
   (onset / offset)
        ↓
Frequency-band detection
        ↓
2D bounding boxes
        ↓
Displayed in GUI
```

The detector first creates a copy of the spectrogram so that preprocessing does not modify the spectrogram displayed by the GUI.

Temporal boundaries are detected first, followed by estimation of the frequency range of each detected syllable.

The resulting detections are converted into the same bounding-box format used for manual annotations.

### Important

Running automatic detection **replaces the currently loaded annotation boxes** with the newly detected boxes.

Therefore, the recommended workflow is:

```text
Load
  ↓
Run Detection
  ↓
Manually inspect detections
  ↓
Add / delete / resize boxes
  ↓
Save
```

If substantial manual annotation has already been performed, save it before running Detection.

---

# Manual Correction

Automatic detection is intended to provide an initial set of candidate USVs.

After detection, the boxes can be manually inspected and corrected using the normal annotation tools.

You can:

* add missed USVs,
* delete false-positive detections,
* move existing boxes,
* resize time boundaries,
* resize frequency boundaries.

Automatic and manually created boxes use the same annotation representation.

---

# Laptop-Optimized Interface

Temp8.1.2 Laptop introduces several layout changes specifically intended for smaller displays.

## Responsive spectrogram

The main spectrogram canvas expands according to the available window space rather than requiring a fixed vertical size.

The time controls also resize horizontally with the window.

## Compact single-USV display

The selected-USV preview uses a compact fixed-height display:

```text
height ≈ 120 px
```

This reduces the amount of vertical space occupied by the preview while retaining access to individual-USV information.

## Frequency axis

The main spectrogram includes a frequency-axis overlay in:

```text
kHz
```

This allows the approximate frequency range of vocalizations and annotation boxes to be read directly from the spectrogram.

## Compact acoustic-feature plots

Instead of placing all acoustic-feature histograms in one horizontal row, the laptop version uses a compact arrangement:

```text
┌──────────┬──────────┬──────────┐
│ Duration │ Loudness │  Purity  │
├──────────┼──────────┼──────────┤
│ Mean Freq│Pitch Var.│          │
└──────────┴──────────┴──────────┘
```

This prevents axis labels and tick labels from overlapping when the GUI is displayed on a narrower screen.

## Air-flow panel removed

The separate air-flow panel used by earlier layouts has been removed.

The available upper display area is instead allocated to the spectrogram, making the vocalization display larger on laptop screens.

---

# Inspecting Individual USVs

Use:

```text
↑ / ↓
```

or the USV index controls to select individual annotations.

The GUI displays:

* Current USV / total USV count
* Individual spectrogram
* Duration
* Loudness
* Spectral purity
* Mean frequency
* Pitch variance

The population distributions of these measurements are displayed in the compact histogram panel.

---

# Saving

Press the **Save** button or:

```text
Ctrl + S
```

The current annotation is saved as:

```text
recording.jld
```

For example:

```text
USV260901_001.h5
        ↓
USV260901_001.jld
```

The bounding-box annotations are stored under:

```julia
"liness"
```

Saving also clears the unsaved-change `*` indicator from the window title.

---

# JPEG Export

Individual USV spectrogram images are automatically exported after saving.

For:

```text
USV260901_001.jld
```

the default output directory is:

```text
USV260901_001_jpegs/
```

The default exported image size is:

```text
128 × 128 px
```

The current exporter uses the low-cut USV spectrogram export workflow with quantile-based contrast adjustment.

These images can subsequently be used for image-based analyses such as clustering, dimensionality reduction, or neural-network models.

---

# Recommended Workflow

For a new recording:

```text
1. Open .h5
       ↓
2. Adjust spectrogram view
       ↓
3. Run Detection
       ↓
4. Inspect detected boxes
       ↓
5. Delete false positives
       ↓
6. Add missed USVs
       ↓
7. Adjust box boundaries
       ↓
8. Save
       ↓
 ┌─────┴─────┐
 ▼           ▼
.jld       JPEGs
```

For an existing annotation:

```text
Open .jld
    ↓
Inspect / edit boxes
    ↓
Save
```

---

# Temp8.1.2 Laptop at a Glance

**Temp8.1.2 Laptop = Temp8.1 detection workflow + laptop-oriented GUI**

The main changes are:

* dedicated `usv_vlaptop.glade` interface,
* responsive spectrogram and GUI controls,
* reduced vertical space usage,
* compact single-USV preview,
* acoustic histograms reorganized into two rows,
* removal of the air-flow panel,
* direct frequency-axis display in kHz,
* improved global keyboard-focus handling.

The automatic detection, manual correction, `.jld` annotation, and JPEG-export workflow are retained.

---

# Notes

Automatic detection should be treated as an **annotation-assistance tool**, not as a substitute for manual inspection.

Detection performance can depend on recording conditions, background noise, microphone configuration, vocalization amplitude, and USV morphology.

For final datasets, visually inspect the automatically generated bounding boxes before saving and downstream analysis.

---

## Current Version

**Temp8.1.2 Laptop**

Main script:

```text
temp8.1.2_laptop.jl
```

GUI definition:

```text
usv_vlaptop.glade
```

Automatic detection functions:

```text
Necessities/usv_detection_corefunctions_temp.jl
```

Other supporting functions are located in:

```text
Necessities/
```
