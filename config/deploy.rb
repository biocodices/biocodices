set :application, 'biocodices'
set :repo_url, 'git@github.com:biocodices/biocodices.git'
set :deploy_to, '~/biocodices/app'
set :keep_releases, 3

set :conda_path, "~/anaconda3/bin"
set :export_path, "export PATH=$PATH:#{fetch(:conda_path)} &&"
set :python, File.join(fetch(:conda_path), "python")

namespace :deploy do
  task :install_app do
    on roles(:app) do
      execute "cd #{current_path} && #{fetch(:python)} setup.py install"
      execute "#{fetch(:export_path)} bioco --version"
    end
  end

  after :deploy, "deploy:install_app"
end

namespace :run do
  task :pipeline, [:base_dir, :options] do |t, args|
    on roles(:app) do
      execute "#{fetch(:export_path)} bioco -d #{args[:base_dir]} #{args[:options]} &"
    end
  end

  task :tail_main_log, [:base_dir] do |t, args|
    on roles(:app) do
      filepath = File.join(args[:base_dir], "pipeline.log")
      execute "tail -n0 -f #{filepath}"
    end
  end

  task :tail_log, [:base_dir, :pattern] do |t, args|
    on roles(:app) do
      pattern1 = File.join(args[:base_dir], "**/*#{args[:pattern]}*.log")
      pattern2 = File.join(args[:base_dir], "*#{args[:pattern]}*.log")
      execute "tail -n0 -f #{pattern1} #{pattern2}"
    end
  end
end
