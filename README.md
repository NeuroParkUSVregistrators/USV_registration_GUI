USV Registration GUI — Temp8.1
A Julia-based graphical interface for manual and automatic registration of mouse ultrasonic vocalizations (USVs) from WaveSurfer .h5 recordings.
Temp8.1 extends the manual USV registration workflow with automatic 2D USV detection, while retaining manual inspection, correction, acoustic feature visualization, annotation saving, and spectrogram-image export.
Features
Load WaveSurfer .h5 recordings
Display USV spectrograms
Automatically detect candidate USVs
Estimate time and frequency boundaries for detected USVs
Convert automatic detections directly into editable GUI bounding boxes
Manually draw additional USV bounding boxes
Move and resize existing annotations
Delete false-positive or unwanted annotations
Save annotations as .jld
Reopen and edit existing .jld annotations
Automatically export individual USVs as JPEG images
Display the current file and unsaved-change status
Navigate recordings using keyboard and mouse controls
Calculate acoustic features for individual and population-level USVs
Displayed acoustic features include:
Duration
Loudness
Spectral purity
Mean frequency
Pitch variance
Repository Structure
USV_registration_GUI/
│
├── temp8.1.jl
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
    ├── usv_detection_corefunctions_temp.jl
    └── JPEGsaver_function_temp.jl
temp8.1.jl is the main GUI script for the automatic-detection-enabled version.
usv_detection_corefunctions_temp.jl contains the core functions used for automatic USV detection.
Requirements
The GUI is written in Julia and currently uses packages including:
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
SparseArrays
Statistics
The current implementation is primarily developed and tested on Windows.
Installation
Clone the repository and switch to the detection-enabled branch:
git clone https://github.com/NeuroParkUSVregistrators/USV_registration_GUI.git
cd USV_registration_GUI
git checkout Detection_activated
Start Julia from the repository directory and activate the environment:
using Pkg
Pkg.activate(".")
Pkg.instantiate()
If a Manifest.toml is provided, Pkg.instantiate() will install the package versions used by the project.
Running the GUI
Run:
include("temp8.1.jl")
The GUI will open a file-selection dialog.
Select either:
recording.h5
for a new registration, or:
recording.jld
to reopen a previously annotated recording.
The corresponding .h5 recording must be available when loading a .jld annotation.
Basic Workflow
WaveSurfer recording (.h5)
          │
          ▼
     Load in GUI
          │
          ▼
  Inspect spectrogram
          │
          ├───────────────┐
          ▼               ▼
 Automatic detection   Manual annotation
          │               │
          └───────┬───────┘
                  ▼
          Inspect / correct
                  │
                  ▼
                Save
               /    \
              ▼      ▼
           .jld     JPEGs
1. Load recording
Open a WaveSurfer .h5 file.
The GUI calculates the spectrogram and displays the recording.
2. Automatic USV detection
Press the Detection button to run automatic USV detection.
The current Temp8.1 detection workflow is:
Spectrogram
    │
    ▼
Denoising
    │
    ▼
Temporal USV detection
    │
    ▼
Frequency-band estimation
    │
    ▼
GUI bounding boxes
Detected USVs are loaded directly into the same annotation structure used for manual registration.
3. Inspect and correct detections
Automatic detections can be edited using the normal GUI controls.
Users can:
delete false-positive detections,
move detected boxes,
resize their time/frequency boundaries, and
manually add missed USVs.
This allows automatic detection to be used as an initial annotation step rather than a replacement for manual inspection.
4. Inspect individual USVs
The selected USV is displayed separately together with its acoustic features.
5. Save
Use the Save button or:
Ctrl + S
The annotation is saved as:
recording.jld
and individual USV spectrogram images are exported automatically.
Automatic Detection
Temp8.1 introduces automatic USV detection integrated directly into the GUI.
The current detection pipeline first denoises a copy of the loaded spectrogram and then identifies temporal onset and offset positions of candidate vocalizations.
For each detected syllable, its local frequency range is subsequently estimated and converted into a two-dimensional GUI bounding box.
Current detection settings include:
Temporal detection
    false-positive mode : confident
    boundary extension  : enabled
    extension SNR       : 1.2
    minimum purity      : 0.15

Frequency detection
    search range        : 15–110 kHz
    selection           : local
    local threshold     : 2.0
    frequency margin    : 4 bins
    minimum box height  : 6 bins
These parameters are currently defined in temp8.1.jl and should be regarded as implementation defaults rather than universal detection thresholds.
Running automatic detection replaces the current annotation list with the detected bounding boxes. Manual correction can then be performed before saving.
Controls
ControlAction
Left mouse drag
Draw new USV
Modified left drag
Move selected USV
Right mouse drag
Resize selected USV
Mouse wheel
Zoom time scale
←
Move backward through recording
→
Move forward through recording
↑
Select next USV
↓
Select previous USV
Ctrl + S
Save
Ctrl + Z
Delete selected USV
Delete
Delete selected USV
Backspace
Delete selected USV
Detection button
Run automatic USV detection
File Status Indicator
The title bar shows the currently loaded filename.
For example:
USV260801_001.h5
When annotations have been changed since the last save:
USV260801_001.h5 *
The * indicates unsaved manual changes.
After saving, the * disappears.
Annotation Files
Annotations are stored in JLD format.
Example:
USV260801_001.jld
Automatic and manual annotations use the same bounding-box representation and can therefore be edited using the same GUI workflow.
The .jld filename corresponds to the original .h5 recording:
USV260801_001.h5
USV260801_001.jld
JPEG Export
Saving also exports individual USV spectrograms.
Default output directory:
USV260801_001_jpegs/
Default image size:
128 × 128
The exported images can subsequently be used for image-based analyses such as dimensionality reduction, clustering, or neural-network models.
Acoustic Measurements
For each registered USV, the GUI calculates:
MeasurementDescription
Duration
Vocalization duration in ms
Loudness
Mean signal power expressed in dB
Spectral purity
Ratio describing concentration of spectral power
Mean frequency
Mean dominant frequency in kHz
Pitch variance
Variance of dominant frequency in log₁₀(Hz²)
Distributions of these measurements are displayed as histograms in the GUI.
Current Version
Temp8.1
Temp8.1 integrates automatic USV detection into the existing manual-registration GUI.
Major features:
Automatic candidate-USV detection
Spectrogram denoising before detection
Temporal onset/offset detection
Local frequency-band estimation for each detected syllable
Automatic conversion of detections into editable 2D GUI bounding boxes
Manual correction of automatic detections
Existing .jld annotation and JPEG-export workflow retained
The intended workflow is therefore:
Automatic detection → Manual review/correction → Save
rather than fully unattended annotation.
Notes
The GUI currently assumes WaveSurfer-style HDF5 recordings and the acquisition structure expected by wavesurfer.jl.
Automatic detection parameters are currently configured for the USV recordings used during development. Detection performance may vary with recording conditions, microphone setup, signal amplitude, background noise, and animal vocalization characteristics.
Automatic detections should therefore be manually inspected before downstream analysis.
Authors / Contributors
NeuroPark USV Registrators
Repository
https://github.com/NeuroParkUSVregistrators/USV_registration_GUI
