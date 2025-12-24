# 100_Trackers_for_1000_Objects

## LEO TLE Dataset – 1000 Objects & 100 Trackers Scenario

### Overview

This dataset contains **1000 Low Earth Orbit (LEO) space objects** represented using **Two-Line Element (TLE)** data, formatted as a CSV file for easy ingestion into analysis and simulation tools.
The dataset is designed to support **Space Situational Awareness (SSA)** studies involving **100 independent trackers monitoring 1000 orbiting objects**.


### Dataset Description
                                                                                                                                   ### Needed more computational power CPU to run this 
* **Number of objects:** 1000
* **Orbit regime:** Low Earth Orbit (LEO)
* **Altitude range:** 450–500 km
* **Orbit type:** Near-circular
* **Inclination range:** ~95°–99° (SSO-like)
* **TLE format:** Standard NORAD Two-Line Element set
* **File format:** CSV


### CSV File Structure

Each row represents one space object.

```csv
csvNORAD_ID,TLE_LINE1,TLE_LINE2
63201,1 63201U 25052P ...,2 63201 97.4217 ...
63202,1 63202U 25053A ...,2 63202 98.1234 ...
```

| Column Name | Description                                        |
| ----------- | -------------------------------------------------- |
| csvNORAD_ID | NORAD catalog identifier extracted from TLE Line 1 |
| TLE_LINE1   | First line of the TLE                              |
| TLE_LINE2   | Second line of the TLE                             |



### 100 Trackers – Concept of Operations (CONOPS)

This dataset is intended to be used in a **100 trackers vs 1000 objects** simulation environment.

**Trackers may include:**

* Ground-based radars
* Ground-based optical sensors
* Space-based optical or IR sensors
* Hybrid SSA sensor networks

**Typical analysis scenarios:**

* Multi-sensor tracking and data association
* Sensor tasking and scheduling
* Coverage and revisit analysis
* Conjunction detection and screening
* Tracking load and scalability testing
* Catalog maintenance simulations

Each tracker can be configured with:

* Field of view constraints
* Pointing laws (nadir, velocity, fixed, target-based)
* Range, elevation, and illumination constraints
* Revisit and detection thresholds



### Supported Tools

This dataset is compatible with:

* **ANSYS STK / ODTK**
* **GMAT**
* **Orekit**
* **Python (sgp4, poliastro, pandas)**
* **MATLAB**
* **Custom SSA simulation frameworks**



### Intended Use Cases

* Space Situational Awareness (SSA) research
* Mission design and analysis
* Sensor performance evaluation
* Constellation traffic modeling
* Academic projects and demonstrations
* Interview assignments and technical assessments



### Notes & Limitations

* TLEs are **synthetically generated** for simulation purposes
* Not intended to represent real, cataloged objects
* Drag, maneuvering, and fragmentation events are not modeled unless explicitly added
* Checksum values are not enforced unless required by the target tool
