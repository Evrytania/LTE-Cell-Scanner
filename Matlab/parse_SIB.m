function parse_SIB(sib_info)

if sib_info.blkcrc ~= 0
    disp('SIB CRC error!');
    return;
end

bits = sib_info.bits;

% 36.331 BCCH-DL-SCH-Message

sib2 = [];
TypeAndInfo.type1 = {sib2};
systemInformation_r8 = {sib-TypeAndInfo};

criticalExtensions.type1 = {systemInformation_r8};
criticalExtensionsFuture = [];
criticalExtensions.type2 = {criticalExtensionsFuture};

systemInformation = {criticalExtensions};

PLMN_Identity = {mcc, mcc};

plmn_IdentityList = 
cellAccessRelatedInfo = {plmn_IdentityList, trackingAreaCode, cellIdentity, cellBarred, intraFreqReselection, csg_Indication, csg_Identity};
systemInformationBlockType1 = {cellAccessRelatedInfo, cellSelectionInfo, p_Max, freqBandIndicator, schedulingInfoList, tdd_Config, si_WindowLength, systemInfoValueTag, nonCriticalExtension};

c1.type1 = {systemInformation};
c1.type2 = {systemInformationBlockType1};

message.type1 = { c1 };
messageClassExtension = [];
message.type2 = { messageClassExtension };

BCCH-DL-SCH-Message = { message };
