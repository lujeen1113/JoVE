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
#include <Python.h>
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
#include "ns3/netanim-module.h"
using namespace std;
using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("QuantumRoutingExample");
vector<vector<double> > readCordinatesFile (std::string node_coordinates_file_name);
int main (int argc, char *argv[])
{
	 std::string node_coordinates_file_name ("src/quantumrouting/examples/4ns.txt");
	vector<vector<double> > coord_array;
	coord_array = readCordinatesFile (node_coordinates_file_name);
	
	
	
	int numNodes = coord_array.size ();
	std::cout<<numNodes<<endl;
        NodeContainer nodeset1;
        nodeset1.Create (numNodes);

        WifiHelper wifi1;
        wifi1.SetStandard (WIFI_PHY_STANDARD_80211b);

        YansWifiPhyHelper wifiPhy1 =  YansWifiPhyHelper::Default ();

        YansWifiChannelHelper wifiChannel1;
        wifiChannel1.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
        wifiChannel1.AddPropagationLoss ("ns3::LogDistancePropagationLossModel");
        wifiPhy1.SetChannel (wifiChannel1.Create ());

        WifiMacHelper wifiMac1;
        wifi1.SetRemoteStationManager ("ns3::ConstantRateWifiManager", "DataMode",StringValue ("DsssRate1Mbps"), "ControlMode",StringValue ("DsssRate1Mbps"));
        wifiMac1.SetType ("ns3::AdhocWifiMac");
        NetDeviceContainer devices1 = wifi1.Install (wifiPhy1, wifiMac1, nodeset1);

        MobilityHelper mobility_nset1;
        Ptr<ListPositionAllocator> positionAlloc_nset1 = CreateObject<ListPositionAllocator> ();

    for (size_t m = 0; m < coord_array.size (); m++)
      {
        positionAlloc_nset1->Add (Vector (coord_array[m][0], coord_array[m][1], 0));
        Ptr<Node> n0 = nodeset1.Get (m);
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
    mobility_nset1.SetPositionAllocator (positionAlloc_nset1);
    mobility_nset1.Install (nodeset1);
    
    
    
  
  QuantumRoutingHelper sr1;
  Ptr<QuantumRoutingTable> srt1 = CreateObject<QuantumRoutingTable> ();
  //Ptr<QuantumRouting> qr= CreateObject<QuantumRouting> ();
  sr1.Set ("RoutingTable", PointerValue (srt1));
  
  InternetStackHelper internet1;
  internet1.SetRoutingHelper (sr1);
  internet1.Install (nodeset1);

  Ipv4AddressHelper ipv4one;
  NS_LOG_INFO ("Assign IP Addresses.");
  ipv4one.SetBase ("10.1.1.0", "255.255.255.0");
  Ipv4InterfaceContainer iface1 = ipv4one.Assign (devices1);
  
  for (int i = 0; i < numNodes; i++)
    {
      nodeset1.Get(i)->setNodeID(i);
	  srt1->AddNode (nodeset1.Get (i), iface1.GetAddress (i));
	  std::cout<<iface1.GetAddress(i)<<endl;
    }
 double flowm[numNodes]={0,1,1,1};
  int trialnum=4;
   double** flowmtrix= new double*[trialnum];
   
   for(int i=0;i<trialnum;i++)
   {
	   flowmtrix[i]=new double[numNodes];
   }
   
   double* flowx=flowm;
   srt1->UpdateRoute (1000000000000);
   srt1->write3pathTable();
   srt1->writeEdgeTable(flowm);
   srt1->getEdgeNum();
   srt1->printStat();
  // srt1->callPythonSolvBundle(true);
   effArray input=srt1->makeEffArray(flowx);
   dwArray output=srt1->makeDWencoding(input);
   
   int srcID=2;
 

   int destID=0;
   
   srt1->callPython2(output);
   srt1->updateRoutingTable();
   
      UdpServerHelper echoServer1 (9);

ApplicationContainer serverApps1 = echoServer1.Install (nodeset1.Get (destID));
serverApps1.Start (Seconds (1.0));
serverApps1.Stop (Seconds (100.0));


UdpEchoClientHelper echoClient1 (iface1.GetAddress (destID), 9);
echoClient1.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient1.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient1.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps1 =
  echoClient1.Install (nodeset1.Get (srcID));
clientApps1.Start (Seconds (2.0));
clientApps1.Stop (Seconds (100.0));

UdpEchoClientHelper echoClient11 (iface1.GetAddress (destID), 9);
echoClient11.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient11.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient11.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps11 =
  echoClient11.Install (nodeset1.Get (srcID+1));
clientApps11.Start (Seconds (2.0));
clientApps11.Stop (Seconds (100.0));

/****************************************************************************************************************************************************/
	
        NodeContainer nodeset2;
        nodeset2.Create (numNodes);

        WifiHelper wifi2;
        wifi2.SetStandard (WIFI_PHY_STANDARD_80211b);

        YansWifiPhyHelper wifiPhy2 =  YansWifiPhyHelper::Default ();

        YansWifiChannelHelper wifiChannel2;
        wifiChannel2.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
        wifiChannel2.AddPropagationLoss ("ns3::LogDistancePropagationLossModel");
        wifiPhy2.SetChannel (wifiChannel2.Create ());

        WifiMacHelper wifiMac2;
        wifi2.SetRemoteStationManager ("ns3::ConstantRateWifiManager", "DataMode",StringValue ("DsssRate1Mbps"), "ControlMode",StringValue ("DsssRate1Mbps"));
        wifiMac2.SetType ("ns3::AdhocWifiMac");
        NetDeviceContainer devices2 = wifi2.Install (wifiPhy2, wifiMac2, nodeset2);

      MobilityHelper mobility_n2;
        Ptr<ListPositionAllocator> positionAlloc_n2 = CreateObject<ListPositionAllocator> ();

    for (size_t m = 0; m < coord_array.size (); m++)
      {
        positionAlloc_n2->Add (Vector (coord_array[m][0], coord_array[m][1], 0));
        Ptr<Node> n0 = nodeset2.Get (m);
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
    mobility_n2.SetPositionAllocator (positionAlloc_n2);
    mobility_n2.Install (nodeset2);
    
   
    
  
  QuantumRoutingHelper sr2;
  Ptr<QuantumRoutingTable> srt2 = CreateObject<QuantumRoutingTable> ();
  //Ptr<QuantumRouting> qr= CreateObject<QuantumRouting> ();
  sr2.Set ("RoutingTable", PointerValue (srt2));
  
  
  InternetStackHelper internet2;
  internet2.SetRoutingHelper (sr2);
  internet2.Install (nodeset2);

  Ipv4AddressHelper ipv4two;
  NS_LOG_INFO ("Assign IP Addresses.");
  ipv4two.SetBase ("10.1.2.0", "255.255.255.0");
  Ipv4InterfaceContainer iface2 = ipv4one.Assign (devices2);
  
  for (int i = 0; i < numNodes; i++)
    {
      nodeset2.Get(i)->setNodeID(i);
	  srt2->AddNode (nodeset2.Get (i), iface2.GetAddress (i));
	  std::cout<<iface2.GetAddress(i)<<endl;
	  std::cout<<nodeset2.Get(i)->GetId();
    }
    

   
  // double* flowx=flowm;
   srt2->UpdateRoute (1000000000000);
   srt2->write3pathTable();
   srt2->writeEdgeTable(flowm);
   srt2->getEdgeNum();
   srt2->printStat();
  // srt2->callPythonSolvBundle(false);
    input=srt2->makeEffArray(flowx);
    output=srt2->makeDWencoding(input);
   
   srt2->callPython2(output);
 
   
   //srt->callPython2(output);
   srt2->updateRoutingTable();
   
      UdpServerHelper echoServer2 (9);

ApplicationContainer serverApps2 = echoServer2.Install (nodeset2.Get (destID));
serverApps2.Start (Seconds (1.0));
serverApps2.Stop (Seconds (100.0));


UdpEchoClientHelper echoClient2 (iface2.GetAddress (destID), 9);
echoClient2.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient2.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient2.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps2 =
  echoClient2.Install (nodeset2.Get (srcID));
clientApps2.Start (Seconds (2.0));
clientApps2.Stop (Seconds (100.0));

UdpEchoClientHelper echoClient21 (iface2.GetAddress (destID), 9);
echoClient21.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient21.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient21.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps21 =
  echoClient21.Install (nodeset2.Get (srcID+1));
clientApps21.Start (Seconds (2.0));
clientApps21.Stop (Seconds (100.0));

/**************************************************************************************************************************************************/

        NodeContainer nodeset3;
        nodeset3.Create (numNodes);

        WifiHelper wifi3;
        wifi3.SetStandard (WIFI_PHY_STANDARD_80211b);

        YansWifiPhyHelper wifiPhy3 =  YansWifiPhyHelper::Default ();

        YansWifiChannelHelper wifiChannel3;
        wifiChannel3.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
        wifiChannel3.AddPropagationLoss ("ns3::LogDistancePropagationLossModel");
        wifiPhy3.SetChannel (wifiChannel3.Create ());

        WifiMacHelper wifiMac3;
        wifi3.SetRemoteStationManager ("ns3::ConstantRateWifiManager", "DataMode",StringValue ("DsssRate1Mbps"), "ControlMode",StringValue ("DsssRate1Mbps"));
        wifiMac3.SetType ("ns3::AdhocWifiMac");
        NetDeviceContainer devices3 = wifi3.Install (wifiPhy3, wifiMac3, nodeset3);

        MobilityHelper mobility_n3;
        Ptr<ListPositionAllocator> positionAlloc_n3 = CreateObject<ListPositionAllocator> ();

    for (size_t m = 0; m < coord_array.size (); m++)
      {
        positionAlloc_n3->Add (Vector (coord_array[m][0], coord_array[m][1], 0));
        Ptr<Node> n0 = nodeset3.Get (m);
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
    mobility_n3.SetPositionAllocator (positionAlloc_n3);
    mobility_n3.Install (nodeset3);
    
    
    
  
  QuantumRoutingHelper sr3;
  Ptr<QuantumRoutingTable> srt3 = CreateObject<QuantumRoutingTable> ();
  //Ptr<QuantumRouting> qr= CreateObject<QuantumRouting> ();
  sr3.Set ("RoutingTable", PointerValue (srt3));
  
InternetStackHelper internet3;
  internet3.SetRoutingHelper (sr3);
  internet3.Install (nodeset3);

  Ipv4AddressHelper ipv4three;
  NS_LOG_INFO ("Assign IP Addresses.");
  ipv4three.SetBase ("10.1.3.0", "255.255.255.0");
  Ipv4InterfaceContainer iface3 = ipv4three.Assign (devices3);
  
  for (int i = 0; i < numNodes; i++)
    {
      nodeset3.Get(i)->setNodeID(i);
	  srt3->AddNode (nodeset3.Get (i), iface3.GetAddress (i));
	//  std::cout<<iface1.GetAddress(i)<<endl;
    }
  
  for (int i = 0; i < numNodes; i++)
    {
      nodeset3.Get(i)->setNodeID(i);
	  srt3->AddNode (nodeset3.Get (i), iface3.GetAddress (i));
	 // std::cout<<iface.GetAddress(i)<<endl;
    }

   
  // double* flowx=flowm;
   srt3->UpdateRoute (1000000000000);
   srt3->write3pathTable();
   srt3->writeEdgeTable(flowm);
   srt3->getEdgeNum();
   srt3->printStat();
  // srt3->callPythonSolvBundle(false);
   
  input=srt3->makeEffArray(flowx);
  output=srt3->makeDWencoding(input);
   
   srt3->callPython2(output);
   
   //srt->callPython2(output);
   srt3->updateRoutingTable();
   
      UdpServerHelper echoServer3 (9);

ApplicationContainer serverApps3 = echoServer3.Install (nodeset3.Get (destID));
serverApps3.Start (Seconds (1.0));
serverApps3.Stop (Seconds (100.0));


UdpEchoClientHelper echoClient3 (iface3.GetAddress (destID), 9);
echoClient3.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient3.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient3.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps3 =
  echoClient3.Install (nodeset2.Get (srcID));
clientApps3.Start (Seconds (2.0));
clientApps3.Stop (Seconds (100.0));

UdpEchoClientHelper echoClient31 (iface3.GetAddress (destID), 9);
echoClient31.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient31.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient31.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps31 =
  echoClient31.Install (nodeset3.Get (srcID+1));
clientApps31.Start (Seconds (2.0));
clientApps31.Stop (Seconds (100.0));


/***************************************************************************************************************************************************/

        NodeContainer nodeset4;
        nodeset4.Create (numNodes);

        WifiHelper wifi4;
        wifi4.SetStandard (WIFI_PHY_STANDARD_80211b);

        YansWifiPhyHelper wifiPhy4 =  YansWifiPhyHelper::Default ();

        YansWifiChannelHelper wifiChannel4;
        wifiChannel4.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
        wifiChannel4.AddPropagationLoss ("ns3::LogDistancePropagationLossModel");
        wifiPhy4.SetChannel (wifiChannel4.Create ());

        WifiMacHelper wifiMac4;
        wifi4.SetRemoteStationManager ("ns3::ConstantRateWifiManager", "DataMode",StringValue ("DsssRate1Mbps"), "ControlMode",StringValue ("DsssRate1Mbps"));
        wifiMac4.SetType ("ns3::AdhocWifiMac");
        NetDeviceContainer devices4 = wifi4.Install (wifiPhy4, wifiMac4, nodeset4);

        MobilityHelper mobility_n4;
        Ptr<ListPositionAllocator> positionAlloc_n4 = CreateObject<ListPositionAllocator> ();

    for (size_t m = 0; m < coord_array.size (); m++)
      {
        positionAlloc_n4->Add (Vector (coord_array[m][0], coord_array[m][1], 0));
        Ptr<Node> n0 = nodeset4.Get (m);
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
    mobility_n4.SetPositionAllocator (positionAlloc_n4);
    mobility_n4.Install (nodeset4);
    
    
    
  
  QuantumRoutingHelper sr4;
  Ptr<QuantumRoutingTable> srt4 = CreateObject<QuantumRoutingTable> ();
  //Ptr<QuantumRouting> qr= CreateObject<QuantumRouting> ();
  sr4.Set ("RoutingTable", PointerValue (srt4));
  
  InternetStackHelper internet4;
  internet4.SetRoutingHelper (sr4);
  internet4.Install (nodeset4);

  Ipv4AddressHelper ipv4four;
  NS_LOG_INFO ("Assign IP Addresses.");
  ipv4four.SetBase ("10.1.4.0", "255.255.255.0");
  Ipv4InterfaceContainer iface4 = ipv4four.Assign (devices4);
  
  for (int i = 0; i < numNodes; i++)
    {
      nodeset4.Get(i)->setNodeID(i);
	  srt4->AddNode (nodeset4.Get (i), iface4.GetAddress (i));
	  std::cout<<iface4.GetAddress(i)<<endl;
    }

   
  // double* flowx=flowm;
   srt4->UpdateRoute (1000000000000);
   srt4->write3pathTable();
   srt4->writeEdgeTable(flowm);
   srt4->getEdgeNum();
   srt4->printStat();
  // srt2->callPythonSolvBundle(false);
   
   input=srt4->makeEffArray(flowx);
   output=srt4->makeDWencoding(input);
   
   srt4->callPython2(output);
   
   //srt->callPython2(output);
   srt4->updateRoutingTable();
   
      UdpServerHelper echoServer4 (9);

ApplicationContainer serverApps4 = echoServer4.Install (nodeset4.Get (destID));
serverApps4.Start (Seconds (1.0));
serverApps4.Stop (Seconds (100.0));


UdpEchoClientHelper echoClient4 (iface4.GetAddress (destID), 9);
echoClient4.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient4.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient4.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps4 =
  echoClient4.Install (nodeset4.Get (srcID));
clientApps4.Start (Seconds (2.0));
clientApps4.Stop (Seconds (100.0));

UdpEchoClientHelper echoClient41 (iface4.GetAddress (destID), 9);
echoClient41.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient41.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient41.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps41 =
  echoClient41.Install (nodeset4.Get (srcID+1));
clientApps41.Start (Seconds (2.0));
clientApps41.Stop (Seconds (100.0));

/*******************************************************************************************************************************************************/

NodeContainer nodeset5;
        nodeset5.Add(nodeset1.Get(0));
        nodeset5.Add(nodeset2.Get(0));
        nodeset5.Add(nodeset3.Get(0));
        nodeset5.Add(nodeset4.Get(0));
        

        WifiHelper wifi5;
        wifi5.SetStandard (WIFI_PHY_STANDARD_80211b);

        YansWifiPhyHelper wifiPhy5 =  YansWifiPhyHelper::Default ();

        YansWifiChannelHelper wifiChannel5;
        wifiChannel5.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel");
        wifiChannel5.AddPropagationLoss ("ns3::LogDistancePropagationLossModel");
        wifiPhy5.SetChannel (wifiChannel5.Create ());

        WifiMacHelper wifiMac5;
        wifi5.SetRemoteStationManager ("ns3::ConstantRateWifiManager", "DataMode",StringValue ("DsssRate1Mbps"), "ControlMode",StringValue ("DsssRate1Mbps"));
        wifiMac5.SetType ("ns3::AdhocWifiMac");
        NetDeviceContainer devices5 = wifi5.Install (wifiPhy5, wifiMac5, nodeset5);

    
  QuantumRoutingHelper sr5;
  Ptr<QuantumRoutingTable> srt5 = CreateObject<QuantumRoutingTable> ();
  //Ptr<QuantumRouting> qr= CreateObject<QuantumRouting> ();
  sr5.Set ("RoutingTable", PointerValue (srt5));
  
  InternetStackHelper internet5;
  internet5.SetRoutingHelper (sr5);
  internet5.Install (nodeset5);

  Ipv4AddressHelper ipv4five;
  NS_LOG_INFO ("Assign IP Addresses.");
  ipv4five.SetBase ("10.1.5.0", "255.255.255.0");
  Ipv4InterfaceContainer iface5 = ipv4five.Assign (devices5);
  
  for (int i = 0; i < numNodes; i++)
    {
      nodeset5.Get(i)->setNodeID(i);
	  srt5->AddNode (nodeset5.Get (i), iface5.GetAddress (i));
	  std::cout<<iface5.GetAddress(i)<<endl;
    }

  // double* flowx=flowm;
   srt5->UpdateRoute (1000000000000);
   srt5->write3pathTable();
   srt5->writeEdgeTable(flowm);
   srt5->getEdgeNum();
   srt5->printStat();
 //  srt5->callPythonSolvBundle(false);
   
    input=srt5->makeEffArray(flowx);
    output=srt5->makeDWencoding(input);
   
   srt5->callPython2(output);
   //srt->callPython2(output);
   srt5->updateRoutingTable();
   
      UdpServerHelper echoServer5 (9);

ApplicationContainer serverApps5 = echoServer5.Install (nodeset5.Get (destID));
serverApps5.Start (Seconds (1.0));
serverApps5.Stop (Seconds (100.0));


UdpEchoClientHelper echoClient5 (iface5.GetAddress (destID), 9);
echoClient5.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient5.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient5.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps5 =
  echoClient5.Install (nodeset5.Get (srcID));
clientApps5.Start (Seconds (2.0));
clientApps5.Stop (Seconds (100.0));
///////////////////////////////////////////////////////////////
UdpEchoClientHelper echoClient51 (iface5.GetAddress (destID), 9);
echoClient51.SetAttribute ("MaxPackets", UintegerValue (1));
echoClient51.SetAttribute ("Interval", TimeValue (Seconds (1.0)));
echoClient51.SetAttribute ("PacketSize", UintegerValue (1024));

ApplicationContainer clientApps51 =
  echoClient51.Install (nodeset5.Get (srcID+1));
clientApps51.Start (Seconds (2.0));
clientApps51.Stop (Seconds (100.0));

  
  AnimationInterface anim("quantumrouting.xml");
  Simulator::Run ();
  Simulator::Destroy ();

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
