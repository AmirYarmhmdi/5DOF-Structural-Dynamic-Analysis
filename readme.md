# 5-Story Building Dynamic Analysis (MATLAB)

This MATLAB project was completed in 2016 as a coursework exercise for the **Structural Dynamics** class.  
It performs a complete **modal analysis** of a 5-story shear building and evaluates its dynamic response under various loading conditions.

## Building Parameters

- Lumped floor mass: 135 kips
- Story stiffness: 315 kips/in
- Target modal damping ratio: 0.03
- Floor heights: [12, 24, 36, 48, 60] ft
- Uniform mass and stiffness distribution across all floors

## Analysis Performed

1. **Modal Analysis**
   - Computation of natural frequencies, periods, and mode shapes
   - Mode shapes plotted for each story under unit amplitude

2. **Step Force Response**
   - Modal response under a step load applied to floors
   - Time history of roof displacement and story drifts
   - Base shear and overturning moment time histories
   - Maximum displacements, drifts, base shear, and overturning moments
   - Graphical visualization of all responses

3. **Harmonic Force Response**
   - Modal response under a sinusoidal force (unit amplitude, first mode frequency)
   - Time history of roof displacement and story drifts
   - Base shear and overturning moment time histories
   - Maximum displacements, drifts, base shear, and overturning moments
   - Graphical visualization of all responses

4. **El-Centro Earthquake Response**
   - Modal response under scaled El-Centro ground acceleration (0.65g)
   - Time history of roof displacement and story drifts
   - Base shear and overturning moment time histories
   - Maximum displacements, drifts, base shear, and overturning moments
   - Both modal superposition and influence/reconstruction methods used
   - Graphical visualization of results

5. **Response Spectrum / Envelope Analysis**
   - Peak responses calculated using modal combination of maximum modal displacements and accelerations
   - Maximum story displacements and drifts
   - Base shear and overturning moment envelopes
   - Graphical visualization

6. **Additional Calculations**
   - Modal effective masses and heights
   - Participation factors for each mode
   - Comparison of maximum floor displacements under El-Centro and envelope methods

## Files

- `5DOF_SDA.m` — Main MATLAB script implementing all calculations
- `elcentro.txt` — El-Centro earthquake acceleration record used in the analysis

## Course Information

- **Course:** Structural Dynamics
- **Professor:** Omid Bahar
- **Student:** Amir Yarmohamadi
- **Year:** 2016
