/* -*-  Mode: C++; c-file-style: "gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2010 University of Arizona
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as 
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Author: Junseok Kim <junseok@email.arizona.edu> <engr.arizona.edu/~junseok>
 */

//    1 -------> 2 -------> 3
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"
#include "ns3/config-store-module.h"
#include "ns3/wifi-module.h"
#include "ns3/internet-module.h"
#include "ns3/quantumrouting-helper.h"
#include "ns3/quantumrouting-table.h"
using namespace std;
using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("QuantumRoutingFQST");
vector<vector<double> > readCordinatesFile (std::string node_coordinates_file_name);
int main (int argc, char *argv[])
{
	 std::string node_coordinates_file_name ("src/quantumrouting/examples/6ns.txt");
	vector<vector<double> > coord_array;
	coord_array = readCordinatesFile (node_coordinates_file_name);
	int numNodes = coord_array.size ();
	std::cout<<numNodes<<endl;
  NodeContainer nodes;
  nodes.Create (numNodes);

  WifiHelper wifi;
  wifi.SetStandard (WIFI_PHY_STANDARD_80211b);

  YansWifiPhyHelper wifiPhy =  YansWifiPhyHelper::Default ();

  YansWifiChannelHelper wifiChannel;
  wifiChannel.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
  wifiChannel.AddPropagationLoss ("ns3::LogDistancePropagationLossModel");
  wifiPhy.SetChannel (wifiChannel.Create ());

  WifiMacHelper wifiMac;
  wifi.SetRemoteStationManager ("ns3::ConstantRateWifiManager", "DataMode",StringValue ("DsssRate1Mbps"), "ControlMode",StringValue ("DsssRate1Mbps"));
  wifiMac.SetType ("ns3::AdhocWifiMac");
  NetDeviceContainer devices = wifi.Install (wifiPhy, wifiMac, nodes);

  MobilityHelper mobility_n;
  Ptr<ListPositionAllocator> positionAlloc_n = CreateObject<ListPositionAllocator> ();

    for (size_t m = 0; m < coord_array.size (); m++)
      {
        positionAlloc_n->Add (Vector (coord_array[m][0], coord_array[m][1], 0));
        Ptr<Node> n0 = nodes.Get (m);
        Ptr<ConstantPositionMobilityModel> nLoc =  n0->GetObject<ConstantPositionMobilityModel> ();
        if (nLoc == 0)
          {
            nLoc = CreateObject<ConstantPositionMobilityModel> ();
            n0->AggregateObject (nLoc);
          }
        // y-coordinates are negated for correct display in NetAnim
        // NetAnim's (0,0) reference coordinates are located on upper left corner
        // by negating the y coordinates, we declare the reference (0,0) coordinate
        // to the bottom left corner
        Vector nVec (coord_array[m][0], -coord_array[m][1], 0);
        nLoc->SetPosition (nVec);

      }
    mobility_n.SetPositionAllocator (positionAlloc_n);
    mobility_n.Install (nodes);
  
  QuantumRoutingHelper sr;
  Ptr<QuantumRoutingTable> srt = CreateObject<QuantumRoutingTable> ();
  sr.Set ("RoutingTable", PointerValue (srt));

  InternetStackHelper internet;
  internet.SetRoutingHelper (sr);
  internet.Install (nodes);

  Ipv4AddressHelper ipv4;
  NS_LOG_INFO ("Assign IP Addresses.");
  ipv4.SetBase ("10.1.1.0", "255.255.255.0");
  Ipv4InterfaceContainer iface = ipv4.Assign (devices);

  for (int i = 0; i < numNodes; i++)
    {
      nodes.Get(i)->setNodeID(i);
	  srt->AddNode (nodes.Get (i), iface.GetAddress (i));
    }
  // double flowx[numNodes]={0,0,0,2,1,0};

   //flowm=flowx;
  /* srt->UpdateRoute (1000000000000);
   srt->write3pathTable();
   srt->writeEdgeTable(flowm);
   srt->getEdgeNum();
   */
   //srt->printStat();
	   int qbuffer[numNodes]={0,0,0,0,0,0};
	  // double flowm[numNodes]={0,0,0,0,0,0};
	   for(int i=1;i<numNodes;i++)
	   {
		   double flowm[numNodes]={0,0,0,0,0,0};
		   flowm[i]=0.2;
		   effArray input=srt->makeEffArray(flowm);
		   dwArray output=srt->makeDWencoding(input);
		   qbuffer[i]=output.qbits;
	   }
	   std::cout<<"qbits against traffic source number"<<endl;
	   for(int i=1;i<numNodes;i++)
	   {
		   std::cout<<qbuffer[i]<<" ";
	   }
	   std::cout<<endl;

   //srt->callPython(output);
   //srt->updateRoutingTable();

  return 0;
}
vector<vector<double> > readCordinatesFile (std::string node_coordinates_file_name)
{
  ifstream node_coordinates_file;
  node_coordinates_file.open (node_coordinates_file_name.c_str (), ios::in);
  if (node_coordinates_file.fail ())
    {
      NS_FATAL_ERROR ("File " << node_coordinates_file_name.c_str () << " not found");
    }
  vector<vector<double> > coord_array;
  int m = 0;

  while (!node_coordinates_file.eof ())
    {
      string line;
      getline (node_coordinates_file, line);

      if (line == "")
        {
          NS_LOG_WARN ("WARNING: Ignoring blank row: " << m);
          break;
        }

      istringstream iss (line);
      double coordinate;
      vector<double> row;
      int n = 0;
      while (iss >> coordinate)
        {
          row.push_back (coordinate);
          n++;
        }

      if (n != 2)
        {
          NS_LOG_ERROR ("ERROR: Number of elements at line#" << m << " is "  << n << " which is not equal to 2 for node coordinates file");
          exit (1);
        }

      else
        {
          coord_array.push_back (row);
        }
      m++;
    }
  node_coordinates_file.close ();
  return coord_array;

}
