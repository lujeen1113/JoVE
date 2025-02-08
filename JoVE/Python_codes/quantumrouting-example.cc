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
  double flowx[numNodes]={0,0,0,2,1,0};
  int trialnum=5;
   double** flowmtrix= new double*[trialnum];
   for(int i=0;i<trialnum;i++)
   {
	   flowmtrix[i]=new double[numNodes];
   }
   
   //double* flowm=flowx;
   srt->UpdateRoute (1000000000000);
   srt->write3pathTable();
   srt->writeEdgeTable(flowx);
   srt->getEdgeNum();
   srt->printStat();
   effArray input=srt->makeEffArray(flowx);
   dwArray output=srt->makeDWencoding(input);
   srt->callPython2(output);
   srt->updateRoutingTable();

  /*effArray input;
  double p[]={1,0,0,1};
  input.penalty=p;
  double inter[][4]={{0,2,0,2.5},{0,0,3,2},{0,0,0,1},{0,0,0,0}};
  double** interaction= new double*[4];
  for(int i=0;i<4;i++)
  {
	  interaction[i]=new double[4];
	  for(int j=0;j<4;j++)
	  {
		  interaction[i][j]=inter[i][j];
	  }
  }
  input.interaction=interaction;
  int var[]={1,3};
  input.var_size=var;
  input.varnum=2;
  input.msize=4;
*/
  //dwArray output=srt->makeDWencoding(input);
//  Py_Initialize();
//  Py_Finalize();
   /*int bitnum=input.msize;
   int flag=1;
   for(int i=0;i<bitnum;i++)
   {
	   if(input.interaction[i][i]!=0)
	   {
		   flag=0;break;
	   }
   }
   std::cout<<"flag is :"<<flag<<endl;

   int num=output.qbits;
   int flag1=1;
   for(int i=0;i<num;i++)
   {
	   if (output.J_core[i][i] !=0)
	   {
		   flag1=0; break;
	   }
   }
   std::cout<<"flag1 :"<<flag1<<endl;
   int flag2=1;
   for(int i=0;i<num;i++)
   {
	   if (output.J_prob[i][i] !=0)
	   {
		   flag2=0; break;
	   }
   }
   std::cout<<"flag2:"<<flag2<<endl;

 std::cout<<"H_core:"<<endl;
  for(int i=0;i<num;i++)
  {
	  std::cout<<output.H_core[i]<<" ";
  }
  std::cout<<endl;

  std::cout<<"H_prob:"<<endl;
   for(int i=0;i<num;i++)
   {
 	  std::cout<<output.H_prob[i]<<" ";
   }
   std::cout<<endl;
   std::cout<<"J_core:"<<endl;
    for(int i=0;i<num;i++)
    {
  	   for(int j=0;j<num;j++){
    	std::cout<<output.J_core[i][j]<<" ";
  	   }
  	   std::cout<<endl;
    }

    std::cout<<"J_prob:"<<endl;
       for(int i=0;i<num;i++)
       {
     	   for(int j=0;j<num;j++){
       	std::cout<<output.J_prob[i][j]<<" ";
     	   }
     	   std::cout<<endl;
       }
*/
  /*double flowm[numNodes]={0,0.3,0,0,0,0};
  srt->UpdateRoute (100);
  srt->write3pathTable();
  srt->writeEdgeTable();
  srt->getEdgeNum();
  srt->printStat();
  effArray str=srt->makeEffArray(flowm);
*/
 /* int bitnum=input.msize;


  std::cout<<"penalty:"<<endl;
   for(int i=0;i<bitnum;i++)
   {
 	  std::cout<<input.penalty[i]<<" ";
   }
   std::cout<<endl;
   std::cout<<"interaction:"<<endl;
    for(int i=0;i<bitnum;i++)
    {
  	   for(int j=0;j<bitnum;j++){
    	std::cout<<input.interaction[i][j]<<" ";
  	   }
  	   std::cout<<endl;
    }
*/
   Ptr<UniformRandomVariable> x = CreateObject<UniformRandomVariable>();
   UdpServerHelper echoServer(9);
   UdpClientHelper clientServer[numNodes];
   ApplicationContainer serverApps1=echoServer.Install (nodes.Get (0));
   serverApps1.Start (Seconds (1.1));
   serverApps1.Stop (Seconds (5.1));


   ApplicationContainer clientApps[numNodes];
   
   int p=0;
   for(double t=0.5;t<5;)
   {
	   double flowx[numNodes]={0,0,0,0,0,0};
	   x->SetAttribute("Min",DoubleValue(0));
	   x->SetAttribute("Max",DoubleValue(numNodes-1));
	   int srcnum=ceil(x->GetValue());
	   //int srcnum=1;
	   int nseqid=0;
	   int options=numNodes-1;
	   double rate=0;
	   for(int i=0;i<srcnum;i++)
	   {   options=options-i;
		   x->SetAttribute("Max",DoubleValue(options));
		   nseqid=ceil(x->GetValue());
		   x->SetAttribute("Max",DoubleValue(5));
		   rate=ceil(x->GetValue());
		   int count=0;
		   int j=1;
		   for(;j<numNodes;j++)
		   {
			   if(flowx[j]==0){count++;}
			   if(count==nseqid){flowx[j]=rate;break;}
		   }

		   //echoServer[i].SetAttribute("Port",UintegerValue(port_num));

		   //clientServer[i](iface.GetAddress (0), 9);
		   clientServer[i].SetAttribute("RemoteAddress",AddressValue(iface.GetAddress(0)));
		   clientServer[i].SetAttribute("RemotePort",UintegerValue(i+numNodes));
		   clientServer[i].SetAttribute ("MaxPackets", UintegerValue (1));
		   clientServer[i].SetAttribute ("Interval", TimeValue (Seconds (0.2)));
		   clientServer[i].SetAttribute ("PacketSize", UintegerValue (1000));
		   clientApps[i]=clientServer[i].Install (nodes.Get (j));
		   clientApps[i].Start (Seconds (0.45+t));
		   clientApps[i].Stop (Seconds (1.43+t));
	   }

	   for(int i=0;i<numNodes;i++)
	   {
           flowmtrix[p][i]=flowx[i];
		   std::cout<<flowmtrix[p][i]<<" ";
	   }
	   std::cout<<endl;
	   EventId Eid=Simulator::Schedule(Seconds(0.44+t),&QuantumRoutingTable::routineCall,srt,flowmtrix,p);
	   std::cout<<Eid.GetUid()<<endl;
	   t=t+1;
	   p=p+1;
   }

  /*UdpServerHelper echoServer1 (9);
  UdpServerHelper echoServer2 (10);
  ApplicationContainer serverApps1 = echoServer2.Install (nodes.Get (0));
  serverApps1.Start (Seconds (1.1));
  serverApps1.Stop (Seconds (5.1));

  UdpClientHelper client1 (iface.GetAddress (0), 10);
  client1.SetAttribute ("MaxPackets", UintegerValue (1));
  client1.SetAttribute ("Interval", TimeValue (Seconds (1)));
  client1.SetAttribute ("PacketSize", UintegerValue (1000));
  ApplicationContainer clientApps;
  clientApps = client1.Install (nodes.Get (3));
  clientApps.Start (Seconds (0.5));
  clientApps.Stop (Seconds (5.0));

  UdpClientHelper client2 (iface.GetAddress(0),10);
  client2.SetAttribute("MaxPackets", UintegerValue(1));
  client2.SetAttribute ("Interval", TimeValue (Seconds (1)));
  client2.SetAttribute ("PacketSize", UintegerValue (1000));
  ApplicationContainer clientApps2;
  clientApps2 =client2.Install(nodes.Get(4));
  clientApps2.Start(Seconds(0.5));
  clientApps2.Stop(Seconds(5.0));*/

  //Ptr<OutputStreamWrapper> routingStream = Create<OutputStreamWrapper> ("quantumrouting.routes",std::ios::out);
  //sr.PrintRoutingTableAllAt(Seconds(5.0),routingStream,Time::S);
  /*UdpServerHelper echoServer2 (10);

  ApplicationContainer serverApps2 = echoServer2.Install (nodes.Get (4));
  serverApps2.Start (Seconds (5.1));
  serverApps2.Stop (Seconds (10.0));

  UdpClientHelper client2 (iface.GetAddress (4), 10);
  client2.SetAttribute ("MaxPackets", UintegerValue (2));
  client2.SetAttribute ("Interval", TimeValue (Seconds (1)));
  client2.SetAttribute ("PacketSize", UintegerValue (1000));


  ApplicationContainer clientApps2;
  clientApps2 = client2.Install (nodes.Get (0));
  clientApps2.Start (Seconds (5.1));
  clientApps2.Stop (Seconds (10.0));
*/

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
