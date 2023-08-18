@echo off

echo %time%

call activate nanodust

(
	start python solo_download_wrap.py "20200615" "20201231"
	start python solo_download_wrap.py "20210101" "20210630"
	start python solo_download_wrap.py "20210701" "20211231"
	start python solo_download_wrap.py "20220101" "20220331"
	start python solo_download_wrap.py "20220401" "20220630"
	start python solo_download_wrap.py "20220701" "20220930"
	start python solo_download_wrap.py "20221001" "20221231"
	start python solo_download_wrap.py "20230101" "20230331"
	start python solo_download_wrap.py "20230401" "20230630"
	start python solo_download_wrap.py "20230701" "20230930"
	
) | set /P "="

echo %time%
