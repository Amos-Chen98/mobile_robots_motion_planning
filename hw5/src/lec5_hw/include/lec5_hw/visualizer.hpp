#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include "lec5_hw/trajectory.hpp"

#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

class Visualizer
{
private:
    ros::NodeHandle nh;

    ros::Publisher wayPointsPub;
    ros::Publisher trajectoryPub;

public:
    Visualizer(ros::NodeHandle &nh_)
        : nh(nh_)
    {
        wayPointsPub = nh.advertise<visualization_msgs::Marker>("/visualizer/waypoints", 10);
        trajectoryPub = nh.advertise<visualization_msgs::Marker>("/visualizer/trajectory", 10);
    }

    template <int D>
    inline void visualize(const Trajectory<D> &traj,
                          const Eigen::Matrix3Xd &route)
    {
        visualization_msgs::Marker wayPointsMarker, trajMarker;

        wayPointsMarker.id = 0;
        wayPointsMarker.type = visualization_msgs::Marker::SPHERE_LIST;
        wayPointsMarker.header.stamp = ros::Time::now();
        wayPointsMarker.header.frame_id = "odom";
        wayPointsMarker.pose.orientation.w = 1.00;
        wayPointsMarker.action = visualization_msgs::Marker::ADD;
        wayPointsMarker.ns = "waypoints";
        wayPointsMarker.color.r = 1.00;
        wayPointsMarker.color.g = 0.00;
        wayPointsMarker.color.b = 0.00;
        wayPointsMarker.color.a = 1.00;
        wayPointsMarker.scale.x = 0.25;
        wayPointsMarker.scale.y = 0.25;
        wayPointsMarker.scale.z = 0.25;

        trajMarker = wayPointsMarker;

        trajMarker.type = visualization_msgs::Marker::LINE_LIST;
        trajMarker.header.frame_id = "odom";
        trajMarker.id = 0;
        trajMarker.ns = "trajectory";
        trajMarker.color.r = 0.00;
        trajMarker.color.g = 0.50;
        trajMarker.color.b = 1.00;
        trajMarker.scale.x = 0.10;

        for (int i = 0; i < route.cols(); ++i)
        {
            geometry_msgs::Point point;
            point.x = route.col(i)(0);
            point.y = route.col(i)(1);
            point.z = route.col(i)(2);
            wayPointsMarker.points.push_back(point);
        }
        wayPointsPub.publish(wayPointsMarker);

        if (traj.getPieceNum() > 0)
        {
            const double T = std::min(0.01, traj.getTotalDuration() / 1000);
            Eigen::Vector3d lastX = traj.getPos(0.0);
            for (double t = T; t < traj.getTotalDuration(); t += T)
            {
                geometry_msgs::Point point;
                Eigen::Vector3d X = traj.getPos(t);
                point.x = lastX(0);
                point.y = lastX(1);
                point.z = lastX(2);
                trajMarker.points.push_back(point);
                point.x = X(0);
                point.y = X(1);
                point.z = X(2);
                trajMarker.points.push_back(point);
                lastX = X;
            }
            trajectoryPub.publish(trajMarker);
        }
        else
        {
            trajMarker.action = visualization_msgs::Marker::DELETEALL;
            trajectoryPub.publish(trajMarker);
        }
    }
};

#endif