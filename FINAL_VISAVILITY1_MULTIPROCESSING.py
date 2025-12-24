"""
FIXED: 100 Space-Based Trackers , 1000 Objects Analysis System
Fixes COM reference errors by proper object lifecycle management

Key Fixes:
1. Create ALL trackers once at start
2. Proper batch cleanup with scenario refresh
3. Smaller batch sizes with memory management
4. COM reference protection
"""

import csv
import tempfile
import os
import datetime as dt
import numpy as np
import gc
import time
from agi.stk12.stkdesktop import STKDesktop
from agi.stk12.stkobjects import *
from agi.stk12.stkutil import *
from agi.stk12.utilities.colors import Color
import json
from multiprocessing import Pool, cpu_count
import multiprocessing as mp

# ============================================================================
# CONFIGURATION
# ============================================================================
CONFIG = {
    'SATELLITES_CSV': 'c:\\Users\\meesa\\Downloads\\stk_100_satellites_orbital_data.csv',
    'OBJECTS_CSV': 'c:\\Users\\meesa\\Downloads\\leo_450_500km_1000_objects_converted.csv',
    'NUM_TRACKERS': 100,
    'NUM_OBJECTS': 1000,
    'ANALYSIS_DAYS': 1,
    'OUTPUT_DIR': 'visibility_results',
    'BATCH_SIZE': 50,  # REDUCED to 10 for better memory management
    'STEP_SIZE': 60,
}

SENSOR_CONFIG = {
    'fov_half_angle': 15.0,
    'max_range': 1000.0,
    'lighting': 'sunlit',
}

import json

def disable_graphics(obj):
    """Completely disable graphics for an object"""
    try:
        obj.Graphics.PassTimes.Visible = False
        obj.Graphics.Show = False
        obj.VO.Visible = False
    except:
        pass

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def ensure_output_dir():
    if not os.path.exists(CONFIG['OUTPUT_DIR']):
        os.makedirs(CONFIG['OUTPUT_DIR'])
        print(f"✓ Created output directory: {CONFIG['OUTPUT_DIR']}")

def load_satellite_trackers(csv_path, num_trackers):
    """Load satellite tracker data with corrected column names"""
    import math
    trackers = []
    try:
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                if i >= num_trackers:
                    break
                
                # Convert True Anomaly to Mean Anomaly
                true_anomaly = float(row.get('TrueAnomaly_deg', 0))
                eccentricity = float(row.get('Eccentricity', 0.0))
                
                e = eccentricity
                TA_rad = math.radians(true_anomaly)
                EA = 2 * math.atan(math.sqrt((1-e)/(1+e)) * math.tan(TA_rad/2))
                MA = EA - e * math.sin(EA)
                mean_anomaly = math.degrees(MA) % 360
                
                trackers.append({
                    'name': row.get('SatelliteName', f'Tracker_{i+1:03d}'),
                    'sma': float(row.get('SMA_km', 6878 + i*10)),
                    'inclination': float(row.get('Inclination_deg', 97.4 + (i % 10) * 0.5)),
                    'eccentricity': eccentricity,
                    'raan': float(row.get('RAAN_deg', (i * 360.0 / num_trackers) % 360)),
                    'arg_perigee': float(row.get('ArgPerigee_deg', (i * 45) % 360)),
                    'mean_anomaly': mean_anomaly
                })
        
        print(f"✓ Loaded {len(trackers)} satellite trackers from {csv_path}")
        return trackers
    except Exception as e:
        print(f"⚠ Could not load satellite CSV: {e}")
        print(f"  Generating {num_trackers} synthetic tracker configurations...")
        return generate_synthetic_trackers(num_trackers)

def generate_synthetic_trackers(num_trackers):
    trackers = []
    for i in range(num_trackers):
        trackers.append({
            'name': f'Tracker_{i+1:03d}',
            'sma': 6878 + (i * 10),
            'inclination': 97.4 + (i % 10) * 0.5,
            'eccentricity': 0.0,
            'raan': (i * 360.0 / num_trackers) % 360,
            'arg_perigee': (i * 45) % 360,
            'mean_anomaly': (i * 30) % 360
        })
    return trackers

def load_tle_objects(csv_path, num_objects):
    tle_data = []
    try:
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            
            for i, row in enumerate(reader):
                if i >= num_objects:
                    break
                
                tle_line1 = row['TLE_LINE1'].strip()
                tle_line2 = row['TLE_LINE2'].strip()
                
                if not tle_line1.startswith('1 ') or not tle_line2.startswith('2 '):
                    print(f"Skipping invalid TLE at row {i+1}")
                    continue
                
                try:
                    norad_id = row['csvNORAD_ID'].strip()
                except:
                    norad_id = tle_line1[2:7].strip()
                
                if len(tle_line1) < 69:
                    tle_line1 = tle_line1.ljust(69)
                if len(tle_line2) < 69:
                    tle_line2 = tle_line2.ljust(69)
                
                tle_data.append({
                    'norad_id': norad_id,
                    'tle_line1': tle_line1,
                    'tle_line2': tle_line2
                })
        
        print(f"✓ Loaded {len(tle_data)} TLE records from {csv_path}")
        
        if tle_data:
            print(f"\nSample TLE:")
            print(f"  NORAD ID: {tle_data[0]['norad_id']}")
            print(f"  Line 1: {tle_data[0]['tle_line1']}")
            print(f"  Line 2: {tle_data[0]['tle_line2']}")
        
        return tle_data
    except Exception as e:
        print(f"ERROR loading TLE data: {e}")
        import traceback
        traceback.print_exc()
        return []

# ============================================================================
# STK OBJECT CREATION
# ============================================================================

def create_tracker_satellite(scenario, config):
    """Create a satellite tracker with sensor - NO GRAPHICS"""
    try:
        tracker = AgSatellite(scenario.Children.New(AgESTKObjectType.eSatellite, config['name']))
        
        # Graphics disabled for speed
        disable_graphics(tracker)

        tracker.SetPropagatorType(AgEVePropagatorType.ePropagatorJ2Perturbation)
        propagator = tracker.Propagator
        keplerian = propagator.InitialState.Representation.ConvertTo(AgEOrbitStateType.eOrbitStateClassical)
        
        keplerian.SizeShapeType = AgEClassicalSizeShape.eSizeShapeSemimajorAxis
        keplerian.SizeShape.SemiMajorAxis = config['sma']
        keplerian.SizeShape.Eccentricity = config['eccentricity']
        
        keplerian.Orientation.Inclination = config['inclination']
        keplerian.Orientation.ArgOfPerigee = config['arg_perigee']
        keplerian.Orientation.AscNodeType = AgEOrientationAscNode.eAscNodeRAAN
        keplerian.Orientation.AscNode.Value = config['raan']
        
        keplerian.LocationType = AgEClassicalLocation.eLocationMeanAnomaly
        keplerian.Location.Value = config['mean_anomaly']
        
        propagator.InitialState.Representation.Assign(keplerian)
        propagator.Propagate()
        
        sensor = AgSensor(tracker.Children.New(AgESTKObjectType.eSensor, f"{config['name']}_Sensor"))
        sensor.CommonTasks.SetPatternSimpleConic(SENSOR_CONFIG['fov_half_angle'], 1.0)
        
        # Graphics disabled for speed
        disable_graphics(sensor)
        
        sensor.SetPointingType(AgESnPointing.eSnPtFixed)
        pointing = sensor.Pointing
        try:
            pointing.Orientation.SetEulerAngles(313, 360, 90, -90)
        except:
            try:
                eulerAngles = pointing.Orientation.ConvertTo(AgEOrientationType.eEulerAngles)
                eulerAngles.A = 360
                eulerAngles.B = 90
                eulerAngles.C = -90
                pointing.Orientation.Assign(eulerAngles)
            except:
                pass
        # add
        constraints = sensor.AccessConstraints
        rangeConstraint = constraints.AddConstraint(AgEAccessConstraints.eCstrRange)
        rangeConstraint.EnableMax = True
        rangeConstraint.Max = SENSOR_CONFIG['max_range']
        
        lightingConstraint = constraints.AddConstraint(AgEAccessConstraints.eCstrLighting)
        lightingConstraint.Condition = 0
        
        return tracker, sensor
    except Exception as e:
        print(f"  ERROR creating tracker {config['name']}: {e}")
        return None, None

def create_space_object(scenario, tle_data):
    """Create a space object from TLE data - NO GRAPHICS"""
    try:
        obj_name = f"Obj_{tle_data['norad_id']}"
        spaceObject = AgSatellite(scenario.Children.New(AgESTKObjectType.eSatellite, obj_name))
        
        # Graphics disabled for speed
        disable_graphics(spaceObject)
        
        spaceObject.SetPropagatorType(AgEVePropagatorType.ePropagatorSGP4)
        sgp4Propagator = spaceObject.Propagator
        
        tle_line1 = tle_data['tle_line1'].strip()
        tle_line2 = tle_data['tle_line2'].strip()
        
        if not tle_line1.startswith('1 ') or not tle_line2.startswith('2 '):
            raise ValueError(f"Invalid TLE format for {obj_name}")
        
        temp_tle = tempfile.NamedTemporaryFile(mode='w', suffix='.tle', delete=False)
        temp_tle.write(tle_line1 + "\n")
        temp_tle.write(tle_line2 + "\n")
        temp_tle.close()
        
        try:
            norad_from_tle = tle_line1[2:7].strip()
            sgp4Propagator.CommonTasks.AddSegsFromFile(norad_from_tle, temp_tle.name)
            sgp4Propagator.Propagate()
        finally:
            os.unlink(temp_tle.name)
        
        return spaceObject
    except Exception as e:
        print(f"  ERROR creating object {tle_data['norad_id']}: {e}")
        try:
            if 'spaceObject' in locals():
                spaceObject.Unload()
        except:
            pass
        return None

def cleanup_objects_by_name(scenario, prefix="Obj_"):
    """Safely clean up objects by name prefix"""
    try:
        objects_to_remove = []
        for i in range(scenario.Children.Count):
            try:
                child = scenario.Children.Item(i)
                if child.InstanceName.startswith(prefix):
                    objects_to_remove.append(child.InstanceName)
            except:
                continue
        
        for obj_name in objects_to_remove:
            try:
                scenario.Children.Unload(AgESTKObjectType.eSatellite, obj_name)
            except:
                pass
        
        return len(objects_to_remove)
    except Exception as e:
        print(f"  Warning during cleanup: {e}")
        return 0

# ============================================================================
# ACCESS COMPUTATION
# ============================================================================

def compute_visibility(sensor, space_objects, scenario, sensor_name):
    """Compute visibility for all objects from one sensor"""
    results = []
    
    for obj in space_objects:
        try:
            access = sensor.GetAccessToObject(obj)
            access.ComputeAccess()
            
            accessDataPrv = access.DataProviders.Item("Access Data").Exec(
                scenario.StartTime, scenario.StopTime
            )
            
            numAccesses = accessDataPrv.Intervals.Count
            
            if numAccesses > 0:
                intervals = []
                for i in range(numAccesses):
                    interval = accessDataPrv.Intervals.Item(i).DataSets
                    intervals.append({
                        'start': interval.GetDataSetByName("Start Time").GetValues()[0],
                        'stop': interval.GetDataSetByName("Stop Time").GetValues()[0],
                        'duration': float(interval.GetDataSetByName("Duration").GetValues()[0])
                    })
                
                results.append({
                    'tracker': sensor_name,
                    'object': obj.InstanceName,
                    'num_passes': numAccesses,
                    'intervals': intervals,
                    'total_duration': sum([a['duration'] for a in intervals])
                })
        except:
            continue
    
    return results


# ============================================================================
# MULTIPROCESSING WORKER
# ============================================================================

def process_single_batch_worker(args):
    """Worker function for multiprocessing - processes one batch"""
    batch_idx, batch_tle, tracker_configs, epoch, end_time, total_batches = args
    
    try:
        # Each worker creates its own STK instance
        uiApp = STKDesktop.StartApplication(visible=False)
        stkRoot = uiApp.Personality2
        try:
            stkRoot.Graphics.Disable()
        except:
            pass
        
        # Create scenario
        stkRoot.NewScenario(f"Batch_{batch_idx}")
        scenario = stkRoot.CurrentScenario
        scenario.StartTime = epoch
        scenario.StopTime = end_time
        stkRoot.Rewind()
        
        # Create trackers
        sensors = []
        for config in tracker_configs:
            tracker, sensor = create_tracker_satellite(scenario, config)
            if tracker and sensor:
                sensors.append({'name': config['name'], 'sensor': sensor})
        
        # Create objects for this batch
        space_objects = []
        for tle in batch_tle:
            obj = create_space_object(scenario, tle)
            if obj:
                space_objects.append(obj)
        
        # Compute visibility
        batch_results = []
        for sensor_info in sensors:
            results = compute_visibility(
                sensor_info['sensor'], 
                space_objects, 
                scenario, 
                sensor_info['name']
            )
            batch_results.extend(results)
        
        # Cleanup
        space_objects.clear()
        gc.collect()
        
        # Close STK instance
        try:
            uiApp.ShutDown()
        except:
            pass
        
        print(f"  ✓ Batch {batch_idx + 1}/{total_batches} complete: {len(batch_results)} detections")
        return batch_results
        
    except Exception as e:
        print(f"  ERROR in batch {batch_idx + 1}: {e}")
        return []


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    start_time = time.time()
    
    print("="*80)
    print("FIXED: 100 SATELLITE TRACKERS × 1000 OBJECTS ANALYSIS")
    print("="*80)
    print(f"\nConfiguration:")
    print(f"  Satellite trackers: {CONFIG['NUM_TRACKERS']}")
    print(f"  Space objects: {CONFIG['NUM_OBJECTS']}")
    print(f"  Analysis period: {CONFIG['ANALYSIS_DAYS']} day(s)")
    print(f"  Batch size: {CONFIG['BATCH_SIZE']} objects/batch")
    
    ensure_output_dir()
    
    # Load data
    print(f"\n{'='*80}")
    print("LOADING DATA")
    print('='*80)
    
    tracker_configs = load_satellite_trackers(CONFIG['SATELLITES_CSV'], CONFIG['NUM_TRACKERS'])
    tle_data = load_tle_objects(CONFIG['OBJECTS_CSV'], CONFIG['NUM_OBJECTS'])
    
    if not tle_data:
        print("ERROR: No TLE data loaded. Exiting.")
        return
    
    print(f"✓ Data loaded: {len(tracker_configs)} trackers, {len(tle_data)} objects")
    
    # Initialize STK
    print(f"\n{'='*80}")
    print("INITIALIZING STK")
    print('='*80)
    try:
        print("Attempting to start STK in background mode...")
        uiApp = STKDesktop.StartApplication(visible=False)
        stkRoot = uiApp.Personality2
        
        # Disable graphics 
        try:
            stkRoot.Graphics.Disable()
        except:
            pass 
        
        print("✓ STK initialized successfully (background mode)")
    except Exception as e:
        print(f"ERROR initializing STK: {e}")
        return
    
    # Create scenario
    print("\nCreating scenario...")
    stkRoot.NewScenario("Multi_Tracker_Visibility")
    scenario = stkRoot.CurrentScenario
    
    epoch = "1 Sep 2025 00:00:00.000"
    end_time = f"{CONFIG['ANALYSIS_DAYS']+1} Sep 2025 00:00:00.000"
    scenario.StartTime = epoch
    scenario.StopTime = end_time
    stkRoot.Rewind()
    
    print(f"✓ Scenario created: {epoch} to {end_time}")
    
   # Process objects in batches using multiprocessing
    print(f"\n{'='*80}")
    print("ANALYZING VISIBILITY IN BATCHES (MULTIPROCESSING)")
    print('='*80)
    
    all_results = []
    total_batches = (len(tle_data) + CONFIG['BATCH_SIZE'] - 1) // CONFIG['BATCH_SIZE']
    
    # Determine number of workers (use 75% of CPU cores)
    num_workers = min(4, max(1, int(cpu_count() * 0.75)))
    print(f"Using {num_workers} parallel workers (CPU cores: {cpu_count()}, limited to 4)")
    
    # Prepare batch arguments
    batch_args = []
    for batch_idx in range(total_batches):
        start_idx = batch_idx * CONFIG['BATCH_SIZE']
        end_idx = min(start_idx + CONFIG['BATCH_SIZE'], len(tle_data))
        batch_tle = tle_data[start_idx:end_idx]
        
        batch_args.append((
            batch_idx,
            batch_tle,
            tracker_configs,
            epoch,
            end_time,
            total_batches
        ))
    
# Process batches in parallel
    print(f"\nProcessing {total_batches} batches in parallel...")
    
    with Pool(processes=num_workers) as pool:
        results_list = pool.map(process_single_batch_worker, batch_args)
    
    # Combine all results
    for batch_results in results_list:
        all_results.extend(batch_results)
    
        print(f"\n✓ All batches complete!")
        print(f"✓ Total detections: {len(all_results)}")
        start_idx = batch_idx * CONFIG['BATCH_SIZE']
        end_idx = min(start_idx + CONFIG['BATCH_SIZE'], len(tle_data))
        batch_tle = tle_data[start_idx:end_idx]
        
        print(f"\nBatch {batch_idx + 1}/{total_batches} (objects {start_idx+1}-{end_idx}):")
        
        # Create space objects for this batch
        print(f"  Creating {len(batch_tle)} space objects...")
        space_objects = []
        for tle in batch_tle:
            obj = create_space_object(scenario, tle)
            if obj:
                space_objects.append(obj)
        
        print(f"  ✓ Created {len(space_objects)} objects")
        
        if len(space_objects) == 0:
            print(f"  ⚠ No objects created in this batch, skipping...")
            continue
        
        # Compute visibility from all trackers
        print(f"  Computing visibility from {CONFIG['NUM_TRACKERS']} trackers...")
        batch_results = []
        
        for i, sensor_info in enumerate():
            if (i + 1) % 20 == 0:
                print(f"    Progress: {i+1}/{CONFIG['NUM_TRACKERS']} trackers processed...")
            
            results = compute_visibility(
                sensor_info['sensor'], 
                space_objects, 
                scenario, 
                sensor_info['name']
            )
            batch_results.extend(results)
        
        all_results.extend(batch_results)
        
        # IMPROVED CLEANUP: Remove objects by name
        print(f"  Cleaning up batch objects...")
        cleaned = cleanup_objects_by_name(scenario, prefix="Obj_")
        print(f"  ✓ Cleaned up {cleaned} objects")
        
        # Clear Python references
        space_objects.clear()
        
        # Force garbage collection
        gc.collect()
        
        # Give STK time to release COM objects
        time.sleep(0.1)
        
        print(f"  ✓ Batch complete. Detections in batch: {len(batch_results)}")
        print(f"  ✓ Total detections so far: {len(all_results)}")
    
    # Analyze results
    print(f"\n{'='*80}")
    print("ANALYZING RESULTS")
    print('='*80)
    
    visible_objects = set([r['object'] for r in all_results])
    
    tracker_stats = {}
    for r in all_results:
        if r['tracker'] not in tracker_stats:
            tracker_stats[r['tracker']] = {'count': 0, 'duration': 0}
        tracker_stats[r['tracker']]['count'] += 1
        tracker_stats[r['tracker']]['duration'] += r['total_duration']
    
    object_stats = {}
    for r in all_results:
        if r['object'] not in object_stats:
            object_stats[r['object']] = {'count': 0, 'duration': 0}
        object_stats[r['object']]['count'] += 1
        object_stats[r['object']]['duration'] += r['total_duration']
    
    # Save results
    print(f"\n{'='*80}")
    print("SAVING RESULTS")
    print('='*80)
    
    output_file = os.path.join(CONFIG['OUTPUT_DIR'], 'visibility_results.json')
    with open(output_file, 'w') as f:
        json.dump({
            'config': CONFIG,
            'num_trackers': CONFIG['NUM_TRACKERS'],
            'num_objects': len(tle_data),
            'num_visible_objects': len(visible_objects),
            'total_detections': len(all_results),
            'tracker_stats': tracker_stats,
            'object_stats': object_stats,
            'all_detections': all_results
        }, f, indent=2)
    
    print(f"✓ Detailed results: {output_file}")
    
    summary_file = os.path.join(CONFIG['OUTPUT_DIR'], 'summary.txt')
    with open(summary_file, 'w') as f:
        f.write("100 SATELLITE TRACKERS × 1000 OBJECTS VISIBILITY ANALYSIS\n")
        f.write("="*80 + "\n\n")
        f.write(f"Analysis Period: {epoch} to {end_time}\n")
        f.write(f"Duration: {CONFIG['ANALYSIS_DAYS']} day(s)\n\n")
        
        f.write("CONFIGURATION:\n")
        f.write(f"  Total trackers: {CONFIG['NUM_TRACKERS']}\n")
        f.write(f"  Total objects analyzed: {len(tle_data)}\n")
        f.write(f"  Sensor FOV: {SENSOR_CONFIG['fov_half_angle']*2}°\n")
        f.write(f"  Max detection range: {SENSOR_CONFIG['max_range']} km\n\n")
        
        f.write("RESULTS:\n")
        if len(visible_objects) > 0:
            f.write(f"  Visible objects: {len(visible_objects)} / {len(tle_data)} ({len(visible_objects)/len(tle_data)*100:.1f}%)\n")
            f.write(f"  Total detections: {len(all_results)}\n")
            f.write(f"  Average detections per visible object: {len(all_results)/len(visible_objects):.1f}\n\n")
        else:
            f.write(f"  Visible objects: 0 / {len(tle_data)} (0.0%)\n")
            f.write(f"  Total detections: 0\n\n")
        
        if tracker_stats:
            f.write("TOP 10 TRACKERS (by detections):\n")
            top_trackers = sorted(tracker_stats.items(), key=lambda x: x[1]['count'], reverse=True)[:10]
            for i, (name, stats) in enumerate(top_trackers, 1):
                f.write(f"  {i}. {name}: {stats['count']} detections, {stats['duration']:.1f}s total\n")
        
        if object_stats:
            f.write("\nTOP 10 OBJECTS (by tracker coverage):\n")
            top_objects = sorted(object_stats.items(), key=lambda x: x[1]['count'], reverse=True)[:10]
            for i, (name, stats) in enumerate(top_objects, 1):
                f.write(f"  {i}. {name}: Seen by {stats['count']} trackers, {stats['duration']:.1f}s total\n")
    
    print(f"✓ Summary: {summary_file}")
    
    elapsed_time = time.time() - start_time
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print('='*80)
    print(f"\nResults Summary:")
    print(f"  ✓ Trackers deployed: {CONFIG['NUM_TRACKERS']}")
    print(f"  ✓ Objects analyzed: {len(tle_data)}")
    if len(visible_objects) > 0:
        print(f"  ✓ Visible objects: {len(visible_objects)} ({len(visible_objects)/len(tle_data)*100:.1f}%)")
        print(f"  ✓ Total detections: {len(all_results)}")
        print(f"  ✓ Avg detections/visible object: {len(all_results)/len(visible_objects):.1f}")
    else:
        print(f"  ✓ Visible objects: 0 (0.0%)")
        print(f"  ✓ Total detections: 0")
    print(f"\n   Total execution time: {elapsed_time:.1f} seconds ({elapsed_time/60:.1f} minutes)")
    print(f"\n   Output directory: {CONFIG['OUTPUT_DIR']}")
    
    print("\n✓ STK running in background mode (no visualization)")
    print("  To visualize results, open the scenario file in STK GUI")  

if __name__ == "__main__":
    main()