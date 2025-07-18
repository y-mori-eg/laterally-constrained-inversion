all: html showcase_plots

distclean: clean
clean:
	rm -f *.png *.gnu *.html *.gif *.webp *.dat epslatex* fit.log soundfit.par *.tmp *.js *.?.svg *.??.svg gnuplot_mouse.css jsdir

GNUPLOT_LIB = ..

#
# Auto-generation of the html demo set using the png terminal.
#
# Several of the demos need manual intervention:
#	fontfile vector epslatex
#
html: airfoil.html animation.html approximate.html armillary.html array.html arrowstyle.html \
      barchart_art.html binary.html binary_polygon.html bins.html bivariat.html boxplot.html \
      boxclusters.html boxes3d.html candlesticks.html chi_shapes.html circles.html cities.html \
      columnhead.html contourfill.html \
      contours.html convex_hull.html custom_contours.html controls.html custom_key.html \
      dashtypes.html datastrings.html dgrid3d.html discrete.html \
      electron.html ellipse.html enhanced_utf8.html epslatex.html errorbars.html \
      fenceplot.html fillbetween.html fillcrvs.html fillstyle.html finance.html fit.html function_block.html \
      heatmaps.html heatmap_4D.html heatmap_points.html \
      hidden.html hidden2.html hidden_compare.html histograms.html histograms2.html hsteps.html \
      histerror.html histogram_colors.html gantt.html image.html image2.html imageNaN.html \
      index.html iris.html iterate.html jitter.html keyentry.html \
      label_stacked_histograms.html layout.html lines_arrows.html linkedaxes.html logic_timing.html \
      map_projection.html margins.html mask_pm3d.html monotonic_spline.html multiaxis.html multiplt.html \
      named_palettes.html nokey.html nonlinear1.html nonlinear2.html nonlinear3.html \
      parallel.html param.html piecewise.html pixmap.html pm3dcolors.html pm3d.html \
      pm3dgamma.html pm3d_clip.html pm3d_lighting.html pointsize.html polar.html polargrid.html \
      poldat.html polar_quadrants.html polygons.html prob2.html prob.html projection.html \
      rainbow.html random.html rank_sequence.html rectangle.html rgba_lines.html \
      argb_hexdata.html rgb_variable.html rotate_labels.html rugplot.html running_avg.html \
      sampling.html scatter.html simple.html singulr.html sectors.html sharpen.html \
      spotlight.html smooth.html smooth_path.html solar_path.html spiderplot.html spline.html \
      smooth_splines.html stats.html steps.html stringvar.html surface1.html surface2.html \
      azimuth.html transparent.html transparent_solids.html textbox.html tics.html \
      timedat.html ttics.html unicode.html using.html \
      varcolor.html vector.html violinplot.html viridis.html voxel.html vplot.html \
      walls.html waterfallplot.html windrose.html world.html zerror.html \
      watchpoints.html isosurface.html

index.html: index.save
	cp index.save index.html
	cp ../../term/js/*.png .

animation.html: animation.html.save
	cp animation.html.save animation.html
	../../src/gnuplot animation.gp

special_functions: index.6new showcase_plots complex_trig.html cerf.html elliptic.html expint.html ibeta.html igamma.html BesselK.html complex_airy.html Dawson.html Fresnel.html lambert.html lnGamma.html synchrotron.html uigamma.html zeta.html
	cp index.6new index.html

circles.html: $(GNUPLOT_LIB)/hidden2.dem ./webify.pl
	export DEMOTERM="pngcairo truecolor font 'arial,8' size 500,300"; ./webify.pl circles

cities.html: $(GNUPLOT_LIB)/cities.dem ./webify.pl
	export DEMOTERM="pngcairo font 'Times,8' fontscale 2.5 size 500,500"; ./webify.pl cities

ellipse.html: random.html

hidden2.html: $(GNUPLOT_LIB)/hidden2.dem ./webify.pl
	export DEMOTERM="pngcairo font 'arial,8' size 540,384"; ./webify.pl hidden2

image.html: $(GNUPLOT_LIB)/image.dem ./webify.pl
	export DEMOTERM="pngcairo font 'arial,8' size 540,384"; ./webify.pl image

image2.html: $(GNUPLOT_LIB)/image2.dem ./webify.pl
	export DEMOTERM="pngcairo font 'arial,8' size 540,384"; ./webify.pl image2

argb_hexdata.html: $(GNUPLOT_LIB)/argb_hexdata.dem ./webify.pl
	export DEMOTERM="pngcairo font 'arial,8' size 540,360"; ./webify.pl argb_hexdata

rgb_variable.html: $(GNUPLOT_LIB)/rgb_variable.dem ./webify.pl
	export DEMOTERM="pngcairo font 'arial,8' size 540,360"; ./webify.pl rgb_variable

rugplot.html: random.html

finance.html: $(GNUPLOT_LIB)/finance.dem ./webify.pl
	export DEMOTERM="pngcairo transparent font 'arial,8' size 660,320"; ./webify.pl finance

linkedaxes.html: $(GNUPLOT_LIB)/linkedaxes.dem ./webify.pl
	export DEMOTERM="pngcairo enhanced font 'arial,9' size 610,320"; ./webify.pl linkedaxes 


multiplt.html: $(GNUPLOT_LIB)/multiplt.dem ./webify.pl
	export DEMOTERM="pngcairo transparent font 'arial,8' size 610,480"; ./webify.pl multiplt

transparent.html: $(GNUPLOT_LIB)/transparent.dem ./webify.pl
	export DEMOTERM="pngcairo transparent truecolor font 'arial,8' size 512,280"; ./webify.pl transparent

transparent_solids.html: $(GNUPLOT_LIB)/transparent_solids.dem ./webify.pl
	export DEMOTERM="pngcairo truecolor font 'arial,9' size 512,384"; ./webify.pl transparent_solids

# The png terminal does not support dashed lines, so let's use pngcairo
lines_arrows.html: $(GNUPLOT_LIB)/lines_arrows.dem ./webify.pl
	export DEMOTERM="pngcairo dash dl 2"; ./webify.pl lines_arrows

rgba_lines.html: $(GNUPLOT_LIB)/rgba_lines.dem ./webify.pl
	export DEMOTERM="pngcairo dash dl 2"; ./webify.pl rgba_lines

complex_trig.html: $(GNUPLOT_LIB)/complex_trig.dem ./webify.pl
	export DEMOTERM="pngcairo enhanced font 'times'"; ./webify.pl complex_trig

showcase_plots: voxel.html pm3d_lighting.html violinplot.html
	cp violinplot.5.png showcase_plot_1.png
	cp voxel.2.png showcase_plot_2.png

%.html: $(GNUPLOT_LIB)/%.dem ./webify.pl
	./webify.pl $*

