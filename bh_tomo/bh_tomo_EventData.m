classdef (ConstructOnLoad) bh_tomo_EventData < event.EventData
   properties
      number = 0;
   end
   methods
      function eventData = bh_tomo_EventData(value)
         eventData.number = value;
      end
   end
end