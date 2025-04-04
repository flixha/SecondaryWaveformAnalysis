
load('arrivals_0.02degGrid_2172events.mat');
arrivals = duplicateArrivalsForDoubleInstrumentSites(arrivals);

%for eclogite transition stuff:
cstation = {'S013','S014','S015','S016','S040','ERM','KRND',...
    'DID','DIDY','DID1','VLY','ATHU','DION','PTL','PSAR','VILL','PROD',...
    'AXAR','EPID','ACOR','LOUT','THAL','EREA','LTK','GUR','KLV','IDHR',...
    'AGG','MAKR','AXAR','ATAL','LKR','DSF','MRKA','THL','KALE',...
    'XOR','NEO','EVR'};

% basically all relevant stations:
cstation = {'ACOR';'AGEO';'AGG';'AGRP';'AIOA';'AKYT';'ALIK';'AMT';'ANDR';'ANKY';...
    'ANX';'AOS';'AT01';'AT02';'AT03';'AT04';'ATAL';'ATH';'ATHU';'AXAR';...
    'AXS';'CHAN';'DID';'DIDY';'DIMT';'DION';'DMLN';'DRO';'DSF';'DSL';...
    'DYR';'EFP';'EPID';'EREA';'ERM';'EVR';'GUR';'HELI';'HORT';'IDHR';'ITM';...
    'KALE';'KARY';'KEAI';'KIMO';'KLMT';'KLV';'KRND';'KTHA';'KTHR';'KYTH';...
    'LAKA';'LIT';'LKR';'LOS';'LOUT';'LTHK';'LTK';'MAKR';'MALA';'MG00';...
    'MHLA';'MHLO';'MRKA';'NEO';'NEST';'NVR';'NYDR';'OHR';'OUR';'PAIG';...
    'PANR';'PDO';'PE01';'PE02';'PE03';'PE04';'PE05';'PE06';'PE07';'PE08';'PE09';...
    'PE10';'PE11';'PENT';'PROD';'PSA1';'PSAM';'PSAR';'PTL';'PVO';'PYL';...
    'PYRG';'RLS';'ROD3';'S001';'S002';'S003';'S004';'S005';'S006';'S007';'S008';...
    'S009';'S010';'S011';'S012';'S013';'S014';'S015';'S016';'S017';'S018';...
    'S019';'S020';'S021';'S022';'S023';'S024';'S025';'S026';'S027';'S028';...
    'S029';'S030';'S031';'S032';'S033';'S034';'S035';'S036';'S037';'S038';'S039';...
    'S040';'S104';'S126';'S129';'SANN';'SANS';'SANT';'SERG';'SERI';'SKIA';...
    'SKY';'SMIA';'SYRO';'TEME';'THAL';'THE';'THL';'TRAZ';'TRI';'TRIP';...
    'TRIZ';'TYRN';'UPR';'VILL';'VLI';'VLMS';'VLS';'VLX';'VLY';'VVK';...
    'XOR';'YDRA'};

newArrivals = arrivals;
for j=1:1:numel(arrivals)
    if isempty(newArrivals(j).arrivals)
        continue
    end
    staFields = fields(newArrivals(j).arrivals);
    for k=1:1:numel(staFields)
        if ~contains(staFields{k},cstation)
            %throw away
            newArrivals(j).arrivals = rmfield(newArrivals(j).arrivals,...
                staFields{k});
        end
    end
end

% % check size of each event's arrival object
% bytesizes = zeros(2172,1);
% for j=1:1:numel(arrivals)
%     A = arrivals(j).arrivals;
%     a = whos('A');
%     bytesizes(j) = a.bytes;
% end

backuparrivals = arrivals;
arrivals=newArrivals;
save('arrivals_0.02degGrid_2172events_relevantStations.mat', 'arrivals');
whos arrivals
whos backuparrivals


% all:
% {'AC01';'AC10';'AC11';'ACOR';'AGEO';'AGG';'AGRP';'AIOA';'AKYT';'ALIK';...
%     'ALNA';'AMT';'ANDR';'ANKY';'ANPA';'ANX';'AOS';'AT01';'AT02';'AT03';...
%     'AT04';'ATAL';'ATH';'ATHU';'AXAR';'AXS';'BIA';'CHAN';'CMBO';'DESF';...
%     'DID';'DIDY';'DIMT';'DION';'DLFA';'DMLN';'DRAG';'DRO';'DSF';'DSL';...
%     'DYR';'EFP';'EPID';'EREA';'ERET';'ERM';'EVGI';'EVR';'FNA';'FOLE';...
%     'FSK';'FYTO';'GRG';'GUR';'HELI';'HORT';'IDHR';'IGT';'IOSI';'IOSP';...
%     'ITHO';'JAN';'KALE';'KARY';'KASA';'KAV';'KAVA';'KBN';'KEAI';'KEF3';...
%     'KEF4';'KEF5';'KEG2';'KEK';'KFL';'KIMO';'KKB';'KLMT';'KLV';'KNS';...
%     'KNT';'KPRO';'KRND';'KRUS';'KTHA';'KTHR';'KTI';'KVR';'KYMI';'KYTH';...
%     'KZN';'LACI';'LAKA';'LIA';'LIT';'LKD';'LKD2';'LKR';'LOS';'LOUT';...
%     'LRSO';'LSK';'LTHK';'LTK';'LXRA';'MAKR';'MALA';'MES2';'MES3';'MET';...
%     'MEV';'MG00';'MGER';'MGNA';'MHLA';'MHLO';'MIL';'MILN';'MILO';'MILS';...
%     'MKIT';'MMB';'MNVA';'MPAR';'MRKA';'MRMA';'MYKO';'N001';'N002';...
%     'N003';'N004';'N005';'N006';'N007';'N008';'N009';'N010';'N011';...
%     'N012';'N013';'N014';'N015';'N016';'N017';'N018';'N019';'N020';...
%     'N021';'N022';'N023';'N024';'N025';'N026';'N027';'N028';'N029';...
%     'N030';'N031';'N033';'N034';'N036';'N037';'N038';'N039';'N040';...
%     'N111';'N136';'N211';'NAIG';'NAXO';'NEAK';'NEO';'NEST';'NSAL';...
%     'NVR';'NYDR';'OHR';'OUR';'PAIG';'PANR';'PARO';'PARS';'PDO';'PE01';...
%     'PE02';'PE03';'PE04';'PE05';'PE07';'PE08';'PE09';'PE10';'PE11';...
%     'PENT';'PHP';'PIL';'PLG';'PROD';'PSA1';'PSAM';'PSAR';'PSDA';'PTL';...
%     'PVO';'PYL';'PYRG';'RGA';'RLS';'ROD3';'RODI';'RODP';'RZN';'S001';...
%     'S002';'S003';'S004';'S005';'S006';'S007';'S008';'S009';'S010';...
%     'S011';'S012';'S013';'S014';'S015';'S016';'S017';'S018';'S019';...
%     'S020';'S021';'S022';'S023';'S024';'S025';'S026';'S027';'S028';...
%     'S029';'S030';'S031';'S032';'S033';'S034';'S035';'S037';'S038';...
%     'S039';'S040';'S104';'S126';'S129';'SANN';'SANS';'SANT';'SAP1';...
%     'SAP2';'SAP3';'SAP4';'SELA';'SERG';'SERI';'SFD';'SGD';'SIFN';'SKIA';...
%     'SKY';'SMIA';'SRN';'STIP';'SYRO';'TEME';'THAL';'THAS';'THE';'THIR';...
%     'THL';'TIR';'TPE';'TRAZ';'TRI';'TRIP';'TRIZ';'TSLK';'TYRN';'UPR';...
%     'VAY';'VIL1';'VIL2';'VILL';'VLI';'VLMS';'VLO';'VLS';'VLX';'VLY';...
%     'VOL2';'VSLK';'VTN';'VVIO';'VVK';'XOR';'XORT';'YDRA';'ZKS'};