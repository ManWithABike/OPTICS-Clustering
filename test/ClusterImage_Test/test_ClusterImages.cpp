// Copyright Ingo Proff 2016.
// https://github.com/CrikeeIP/OPTICS-Clustering
// Distributed under the MIT Software License (X11 license).
// (See accompanying file LICENSE)

#include "../../include/optics/optics.hpp"
#include "../../include/optics/bgr_image.hpp"

#include <vector>




inline std::vector<std::array<int, 2>> load_points_from_image( const std::string& file_path )
{
	const bgr_image img = imread( file_path );
	assert( img.size().area() != 0 );
	std::vector<std::array<int, 2>> result;
	for ( std::size_t y = 0; y < img.size().height_; ++y )
	{
		for ( std::size_t x = 0; x < img.size().width_; ++x )
		{
			if ( img.pix( img_pos(x, y) ) != bgr_col(255,255,255) )
			{
				result.push_back( { static_cast<int>(x), static_cast<int>(y) } );
			}
		}
	}
	return result;
}


int main__()
{
	std::string path = "./test/ClusterImage_Test/ClusterImage_1.ppm";
	const auto points = load_points_from_image( path );
	std::cout << "Extracted " << points.size() << " points from the image" << std::endl;

	std::size_t min_pts = 15;
	auto eps = 2 * optics::epsilon_estimation( points, min_pts );
	std::cout << "Estimated epsilon = " << eps << std::endl;

	std::cout << "Computing Reachability-Distances..." << std::endl;
	const auto reach_dists = optics::compute_reachability_dists<1000>( points, min_pts, eps );
	std::cout << "Done!" << std::endl;

	//std::string str = fplus::show(reach_dists);
	//std::cout << "Written: " << fplus::write_text_file("ReachDists_Text.txt", str)() << std::endl;

	/*for ( std::size_t l = 0; l < 1000; l++ ) {
		auto reach_dists_ = optics::compute_reachability_dists( points, min_pts, eps );
		if ( l % 10 == 0 ) std::cout << l << std::endl;
	}*/

	std::cout << "Drawing reachability-plot" << std::endl;
	auto reach_img = optics::draw_reachability_plot( reach_dists );
	reach_img.save( "./ReachDists" );

   std::cout << "Drawing 2D-Clusters" << std::endl;
	{
		auto clusters = optics::get_cluster_points( reach_dists, 10, points );
		auto img = optics::draw_2d_clusters( clusters );
		img.save( "./ClusterImg_10" );
	}
	{
		auto clusters = optics::get_cluster_points( reach_dists, 20, points );
		auto img = optics::draw_2d_clusters( clusters );
		img.save( "./ClusterImg_20" );
	}

   std::cout << "Computing Chi-Clusters" << std::endl;
	{
		//Chi Clusters
		double chi = 0.02;
		double steep_area_min_diff = 0.25;

		auto chi_img = optics::draw_reachability_plot_with_chi_clusters( reach_dists, chi, min_pts, steep_area_min_diff );
		chi_img.save( "./Reachdists_ChiCluster" );

		auto chi_clusters = optics::get_chi_clusters( reach_dists, chi, min_pts, steep_area_min_diff );
		std::vector<optics::chi_cluster_indices> chi_clusters_flat;
		chi_clusters_flat = fplus::concat( fplus::transform( optics::flatten_dfs<optics::chi_cluster_indices>, chi_clusters ) );
		auto cluster_pts = optics::get_cluster_points( reach_dists, chi_clusters_flat, points );
		auto chi_cluster_img = optics::draw_2d_clusters( cluster_pts );
		chi_cluster_img.save( "./ChiClusterImg" );
	}

	return 0;
}
