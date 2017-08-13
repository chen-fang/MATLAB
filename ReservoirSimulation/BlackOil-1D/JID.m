function [f] = JID( PhaseID, CellID ) % 2-phase

f = 2*( CellID - 1 ) + PhaseID;