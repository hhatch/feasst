***************************
List of all capabilities
***************************

FEASST simulations are conducted by creating a series of classes in a specific order.
The first word in each line of a FEASST input text file is the name of a class, followed by pairs of class arguments.

The basic classes are as follows.
:cpp:class:`Model <feasst::Model>` computes the potential energy between two sites.
To find all available :cpp:class:`Models <feasst::Model>`, look at the subclasses in the documentation for :cpp:class:`one <feasst::ModelOneBody>` or :cpp:class:`two <feasst::ModelTwoBody>` body Models.
:cpp:class:`VisitModel <feasst::VisitModel>` loops a :cpp:class:`Model <feasst::Model>` over all selected particles and sites, or contains a many-body potential such as :cpp:class:`LongRangeCorrections <feasst::LongRangeCorrections>` or :cpp:class:`Ewald <feasst::Ewald>`.
:cpp:class:`Potential <feasst::Potential>` contains both :cpp:class:`Model <feasst::Model>` and :cpp:class:`VisitModel <feasst::VisitModel>`.
:cpp:class:`Trial <feasst::Trial>` is a :cpp:class:`MonteCarlo <feasst::MonteCarlo>` trial move.
:cpp:class:`Criteria <feasst::Criteria>` determine if the :cpp:class:`Trial <feasst::Trial>` is accepted or rejected.
:cpp:class:`FlatHistogram <feasst::FlatHistogram>` methods :cpp:class:`Bias <feasst::Bias>` along a :cpp:class:`Macrostate <feasst::Macrostate>`.
:cpp:class:`Analyze <feasst::Analyze>` and :cpp:class:`Modify <feasst::Modify>` occur every fixed number of :cpp:class:`Trials <feasst::Trial>`, but an :cpp:class:`Action <feasst::Action>` happens only once.
These base classes have all their subclasses listed in their respective documentation.

A FEASST plugin is a collection of related classes.
Unnecessary plugins can be removed from the FEASST_PLUGINS variable in /path/to/feasst/CMakeLists.txt.
Similarly, plugins can also be added this way.
In both cases, FEASST must be reinstalled.

See the example plugin as a template for creating your own class or plugin.

List of plugins
--------------------------------

The plugins and classes listed below represent all the publicly available capabilities of FEASST.

.. toctree::

   threads/README
   utils/README
   math/README
   configuration/README
   system/README
   monte_carlo/README
   models/README
   steppers/README
   flat_histogram/README
   patch/README
   mayer/README
   xtc/README
   chain/README
   shape/README
   confinement/README
   charge/README
   opt_lj/README
   cluster/README
   egce/README
   morph/README
   beta_expanded/README
   prefetch/README
   aniso/README
   example/README
   fftw/README
   netcdf/README

