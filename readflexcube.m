function [ cube, meta ] = readflexcube(fname, varargin)
% READFLEXCUBE reads XYZC data from a PerkinElmer Flex file
% [ CUBE, META ] = READFLEXCUBE(FNAME)
% [ CUBE, META ] = READFLEXCUBE(FNAME, 'sizeZ' sizeZ, 'sizeC', sizeC)
% [ ~, META ] = READFLEXCUBE(FNAME, 'ImageData', false)
%
%   FNAME is the filename of the flex file
%
%   CUBE is a 4-D array containing image data
%   META is a struct with image metadata
%
% Two optional keyword arguments can be used to overwrite the
% auto-detection of the stack geometry:
%
%   sizeZ is the number of Z stacks in the flex file
%   sizeC is the number of channels
%
% The keyword 'ImageData' can be used to control whether image data should
% be read (otherwise only metadata is returned). Defaults to true.
%
% Author: Bob Pepin (bobpepin@gmail.com)

parser = inputParser;
addParamValue(parser, 'Info', []);
addParamValue(parser, 'sizeZ', []);
addParamValue(parser, 'sizeC', []);
addParamValue(parser, 'PixelSizeZ', []);
addParamValue(parser, 'ImageData', true);
parse(parser, varargin{:});
args = parser.Results;

[~, f] = fileattrib(fname);
fname = f.Name;

info = [];
if nargout > 1
    [meta, info] = flexfileinfo(fname);
     
    if not(isempty(args.PixelSizeZ))
        pixelZ = args.PixelSizeZ;
    else
        pixelZ = nan;
    end
    
    meta.pixelSizeZ = pixelZ;
end

if ~args.ImageData
    cube = [];
    return;
end

if not(isempty(args.Info))
    info = args.Info;
elseif isempty(info)
    info = imfinfo(fname);
end

sizeX = info(1).Width;
sizeY = info(1).Height;
sizeZ = 1;
sizeC = 1;

if not(isempty(args.sizeZ))
    sizeZ = args.sizeZ;
end

if not(isempty(args.sizeC))
    sizeC = args.sizeC;
end

totalcount = numel(info);

if totalcount ~= sizeZ*sizeC
    if isempty(args.sizeZ)
       sizeZ = floor(totalcount / sizeC);
    elseif isempty(args.sizeC)
       sizeC = floor(totalcount / sizeZ);
    end
end

cube = zeros([sizeY(), sizeX(), sizeZ, sizeC], 'uint16');

% getZCTCoords = @(i) [mod(i, sizeZ), ...
%                      mod(floor(i / sizeZ), sizeC), ...
%                      floor(floor(i / sizeZ) / sizeC)];

getCZTCoords = @(i) [mod(floor(i / sizeC), sizeZ), ...
                     mod(i, sizeC), ...
                     floor(floor(i / sizeC) / sizeZ)];

w = warning('off','all');
tf = Tiff(fname);
for i=1:totalcount
   zct = getCZTCoords(i - 1) + 1;
   cube(:, :, zct(1), zct(2)) = tf.read();
   if i < totalcount
       tf.nextDirectory();
   end
end
warning(w);



end

function [flexinfo, info] = flexfileinfo(fname)

flexinfo.filename = fname;

retries = 5;
while retries > 0, try %#ok<ALIGN>
    info = imfinfo(fname);
    retries = 0;
catch err
    retries = retries - 1;
    if retries == 0
        rethrow(err)
    end
    fprintf('Unable to open file %s, %d retries left...\n', fname, retries);
end, end

xmlstr = info(1).UnknownTags(1).Value;
xml = xmlreadstring(xmlstr);

flexinfo.areaname = char(xml.getElementsByTagName('AreaName').item(0).getTextContent());
attr = xml.getElementsByTagName('WellCoordinate').item(0).getAttributes();
flexinfo.row = str2double(char(attr.getNamedItem('Row').getNodeValue()));
flexinfo.column = str2double(char(attr.getNamedItem('Col').getNodeValue()));
flexinfo.meas = str2double(char(xml.getElementsByTagName('OperaMeasNo').item(0).getTextContent()));
flexinfo.field = str2double(char(xml.getElementsByTagName('Images').item(0).getElementsByTagName('Sublayout').item(0).getTextContent()));
flexinfo.timestamp = char(xml.getElementsByTagName('DateTime').item(0).getTextContent());

flexinfo.pixelSizeX = str2double(char(xml.getElementsByTagName('ImageResolutionX').item(0).getTextContent()));
flexinfo.pixelSizeY = str2double(char(xml.getElementsByTagName('ImageResolutionY').item(0).getTextContent()));

flexinfo.barcode = char(xml.getElementsByTagName('Barcode').item(0).getTextContent());
flexinfo.plateColumns = str2double(char(xml.getElementsByTagName('YSize').item(0).getTextContent()));
flexinfo.plateRows = str2double(char(xml.getElementsByTagName('XSize').item(0).getTextContent()));

end

function [parseResult,p] = xmlreadstring(stringToParse,varargin)
%XMLREADSTRING Modified XMLREAD function to read XML data from a string.
% Author: Luis Cantero.
% The MathWorks.

p = locGetParser(varargin);
locSetEntityResolver(p,varargin);
locSetErrorHandler(p,varargin);

% Parse and return.
parseStringBuffer = java.io.StringBufferInputStream(stringToParse);
parseResult = p.parse(parseStringBuffer);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = locGetParser(args)

p = [];
for i=1:length(args)
    if isa(args{i},'javax.xml.parsers.DocumentBuilderFactory')
        javaMethod('setValidating',args{i},locIsValidating(args));
        p = javaMethod('newDocumentBuilder',args{i});
        break;
    elseif isa(args{i},'javax.xml.parsers.DocumentBuilder')
        p = args{i};
        break;
    end
end

if isempty(p)
    parserFactory = javaMethod('newInstance',...
        'javax.xml.parsers.DocumentBuilderFactory');
        
    javaMethod('setValidating',parserFactory,locIsValidating(args));
    %javaMethod('setIgnoringElementContentWhitespace',parserFactory,1);
    %ignorable whitespace requires a validating parser and a content model
    p = javaMethod('newDocumentBuilder',parserFactory);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf=locIsValidating(args)

tf=any(strcmp(args,'-validating'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function locSetEntityResolver(p,args)

for i=1:length(args)
    if isa(args{i},'org.xml.sax.EntityResolver')
        p.setEntityResolver(args{i});
        break;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function locSetErrorHandler(p,args)

for i=1:length(args)
    if isa(args{i},'org.xml.sax.ErrorHandler')
        p.setErrorHandler(args{i});
        break;
    end
end

end

