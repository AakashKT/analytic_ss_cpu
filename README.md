Fast Analytic Soft Shadows from Area Lights, EGSR 2021
===============
<b>Aakash KT <sup><a href="http://cvit.iiit.ac.in/">[1]</a></sup></b>, Parikshit Sakurikar <sup><a href="https://dreamvu.com/">[1, 2]</a></sup>, P. J. Narayanan <sup><a href="http://cvit.iiit.ac.in/">[1]</a></sup>
<br>
<span>
    <a target="_blank" href="#">[Publisher's Version]</a>
    <a target="_blank" href="https://drive.google.com/file/d/1Ssgnxu-yjYeDAasBKEuaM6i9q6KcK5g1/view">[Paper]</a>
    <a target="_blank" href="https://aakashkt.github.io/analytic_ss_bibtex.txt">[bibtex]</a>
</span>
<br><br>
<b>Abstract</b>
<p>
In this paper, we present the first method to analytically compute shading and soft shadows for physically based BRDFs from
arbitrary area lights. We observe that for a given shading point, shadowed radiance can be computed by analytically integrating over the visible portion of the light source using Linearly Transformed Cosines (LTCs). We present a structured approach to project, re-order and horizon-clip spherical polygons of arbitrary lights and occluders. The visible portion is then computed by multiple repetitive set difference operations. Our method produces noise-free shading and soft-shadows and outperforms ray-tracing within the same compute budget. We further optimize our algorithm for convex light and occluder meshes by projecting the silhouette edges as viewed from the shading point to a spherical polygon, and performing one set difference operation thereby achieving a speedup of more than 2×. We analyze the run-time performance of our method and show rendering results on several scenes with multiple light sources and complex occluders. We demonstrate superior results compared to prior work that uses analytic shading with stochastic shadow computation for area lights.
</p>

Code
---------
This code is written on top of PBRT-v3. Our analytic method is implemented as integrator plugins. The main code resides in the following files:<br>
<ul>
<li><b>src/integrators/ltc.cpp:</b> LTC implementation</li>
<li><b>src/integrators/ratio.cpp:</b> Ratio Estimator (See <a href="https://research.nvidia.com/publication/2018-05_Combining-Analytic-Direct">here</a>.)</li>
<li><b>src/integrators/ltc+shadow.cpp:</b> Our method (Per-triangle)</li>
<li><b>src/integrators/ltc+silhouette+shadow.cpp:</b> Our method (Silhouette edges)</li>
<li><b>src/integrators/polygon*.cpp:</b> Set difference using Greiner-Hormann, adapted from <a href="https://github.com/Lecanyu/PolygonClipping">here</a></li>
</ul>

No extra dependencies are required, hence you can build this code following the original instructions from <a href="https://github.com/mmp/pbrt-v3">PBRT-v3</a>.
<br>
You can find scenes used in the paper <a href="https://iiitaphyd-my.sharepoint.com/:f:/g/personal/aakash_kt_research_iiit_ac_in/Emm-eVr3AodEpjMEr8XPMCIBdxOVABKX2JdhbWDpRwpGRA?e=zlxVTd">here</a>.

Rendering & Visualizing
---
To render using our silhouette edge method, first change the sampler and integrator in your scene file like so:
```
Sampler "stratified" "integer xsamples" 1 "integer ysamples" 1 "bool jitter" "false"
Integrator "ltc+silhouette+shadow"
```
The various integrator names correspond exactly to the file names given above. Now run PBRT as:
```
./pbrt scene.pbrt -outfile render.png
```

To visualize spherical polygons and difference polygons, first uncomment all lines between "LOG" and "ENDLOG" in the file of the integrator of your choice and rebuild PBRT. Then run PBRT with one thread:
```
./pbrt scene.pbrt -outfile render.png --nthreads 1
```
This will save a "pixelOutput.txt" in the build directory. To visualize the polygons, use the "visualize.py" script like so:
```
python visualize.py --image_file render.png --pixel_file pixelOutput.txt
```
This will give you an interactive polygon visualizer, where you can click on any point in the image and all polygons (blue for light, red for occluder) will be visualized in the middle and difference polygons will be visualized in the right. Press "B" to toggle between hiding and showing occluders, and press "H" to toggle between a lat-long (spherical) projection and a 2D (planar) projection.

To render with the ratio estimator, prepare three files corresponding to Sn, Un, and U (see provided scenes, scene_un.pbrt, scene_sn.pbrt, scene_ltc.pbrt respectively). Render with ratio estimator like so:
```
python analytic_stochastic.py --scene_dir /path/to/scene/files
```
This will render and save individual three outputs and the composite output, using bilateral filter for denoising. All timings will also be reported.

Acknowledgements
------------
We thank the reviewers for their valuable feedback and insights. This research was partially funded by the "Kohli Fellowship" of KCIS. <br>
We also thank Matt Pharr, Wenzel Jakob, and Greg Humphreys for their amazing book <i>Physically Based Rendering: From Theory to Implementation</i>.


