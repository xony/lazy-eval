ACLOCAL_AMFLAGS = -I m4


AM_YFLAGS=-d
AM_CXXFLAGS = -fPIC -I../sequoia-core/src -I../sequoia-core/contrib/sparsehash-2.0.2/src
AM_CPPFLAGS =                   $(BOOST_CPPFLAGS)
#CXXFLAGS += -DNDEBUG
#CXXFLAGS += -O0 -DDODEBUG -pg -g 
AM_LDFLAGS = -lfl -ldl -lpthread $(BOOST_SYSTEM_LDFLAGS)
LDADD  = ../sequoia-core/src/libsequoia.la \
	$(BOOST_SYSTEM_LIBS) \
                                $(BOOST_GRAPH_LIBS) \
                                $(BOOST_TIMER_LIBS) \
                                $(BOOST_CHRONO_LIBS) \
                                $(BOOST_FILESYSTEM_LIBS) \
                                $(BOOST_THREAD_LIBS)

		

bin_PROGRAMS =  gen_clique_tdc \
		vertex_cover_lazy \
		vertex_cover_lazy_refactored  \
		sat_to_vc \
		vertex_cover \
		vertex_cover_lazyonly \
		gen_tdc \
		vertex_cover_lazy_mindegree
		

		 


vertex_cover_SOURCES = vertex_cover.cpp
vertex_cover_lazyonly_SOURCES = vertex_cover_lazyonly.cpp
vertex_cover_lazyonly_LDFLAGS  = $(AM_LDFLAGS) -lgmp -lgmpxx
vertex_cover_lazy_SOURCES = vertex_cover_lazy.cpp
vertex_cover_lazy_LDFLAGS  = $(AM_LDFLAGS) -lgmp -lgmpxx
vertex_cover_lazy_refactored_SOURCES = vertex_cover_lazy_refactored.cpp
vertex_cover_lazy_refactored_LDFLAGS  = $(AM_LDFLAGS) -lgmp -lgmpxx
vertex_cover_lazy_mindegree_SOURCES = vertex_cover_lazy_mindegree.cpp
vertex_cover_lazy_mindegree_LDFLAGS  = $(AM_LDFLAGS) -lgmp -lgmpxx
#gen_grid_tdc_SOURCES =  gen_grid_tdc.cpp
#gen_rand_graph_SOURCES =  gen_rand_graph.cpp
gen_clique_tdc_SOURCES =  gen_clique_tdc.cpp
sat_to_vc_SOURCES =  sat_to_vc.cpp
gen_tdc_SOURCES =  gen_tdc.cpp
