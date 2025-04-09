🐛 Termite Tandem Behavior Analysis
This project investigates tandem running behavior in termites, using high-resolution position data extracted from video tracking. The analysis focuses on identifying tandem runs, measuring their durations, detecting leadership (male or female), and analyzing how behaviors vary by arena size.

📁 Project Structure
bash
Copy
Edit
project_root/
│
├── analysis_WC250331.R       # Main analysis pipeline
├── data_fmt/
│   ├── df_all.rda            # Full behavior-enriched dataset
│   ├── df_body_scaled.rda    # Scaled body size data (mm)
│   ├── df_lead.rda           # Processed leadership behavior
│   └── data_raw_df.feather   # Raw tracking data (from SLEAP)
├── data_raw_bodysize.csv     # Raw body measurements
└── README.md                 # This file
🔧 Setup & Dependencies
This project is written in R. You’ll need the following packages:

r
Copy
Edit
install.packages(c("stringr", "data.table", "arrow", "dplyr", "MASS", "ggplot2",
                   "patchwork", "knitr", "survival", "survminer", "zoo",
                   "cowplot", "coxme", "tidyr"))
🧠 Research Goals
Determine how sex (male vs. female) affects tandem leadership.

Evaluate how arena size (90mm vs. 150mm) influences behavior.

Understand duration and stability of tandem runs and separation events.

Quantify switches in leadership during tandem.

📊 Workflow Summary
1. Data Processing
Scale raw tracking data (SLEAP) from pixels to millimeters using arena-specific scaling.

Downsample frames to 5 FPS for analysis consistency.

Merge position data with body size metrics.

2. Tandem Identification
A tandem run is defined by:

Speed threshold

Distance between individuals

Angular alignment

Tandem smoothing removes short spurious events.

3. Leader Detection
Leaders are identified based on positional geometry:

Whichever termite's head is closer to the partner's abdomen is the follower.

Smoothing is applied to leadership frames for reliability.

4. Feature Engineering
Calculate durations of:

Tandem runs (split by sex and dish size)

Separation events (post-tandem)

Count switches in leadership.

Label who was leading prior to separation.

5. Statistical Modeling
Mixed-effects Cox proportional hazards models (coxme) analyze:

Whether male or female leaders have longer tandem durations.

Differences in behavior between dish sizes.

6. Visualization
Survival analysis plots:

Tandem duration by sex

Separation duration by previous leader

Comparisons by arena size

Boxplots:

Speed distributions by sex and dish size

Leadership switching rates

📌 Key Parameters
r
Copy
Edit
tandemAngle <- 60 * (pi / 180)  # Angular difference threshold
tandemsmooth <- 20              # Smoothing window for tandem detection
leadsmooth <- 15                # Smoothing window for leadership
tandemSpeed <- 1.213            # Minimum speed (mm/s) to count as tandem
🐜 Data Sources
Tracking Data: data_raw_df.feather from SLEAP.

Termite Body Sizes: data_raw_bodysize.csv

Processed Data: Includes .rda outputs with leadership and behavior metadata.

✅ Status
 Preprocessing complete

 Tandem behavior annotated

 Statistical modeling done

 Survival visualizations created

 Final publication-ready plots and write-up in progress
